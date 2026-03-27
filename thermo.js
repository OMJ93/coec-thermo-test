async function loadThermoDataset(jsonPath = "./thermo_dataset.json") {
  const response = await fetch(jsonPath);
  if (!response.ok) {
    throw new Error(`Failed to load dataset: ${response.status}`);
  }
  return await response.json();
}

function nearlyEqual(a, b, tol = 1e-9) {
  return Math.abs(Number(a) - Number(b)) <= tol;
}

function filterRowsByCondition(dataset, {
  T_C,
  ratio_h2o_to_co2,
  conversion,
  air_outlet_o2_n2_ratio,
}) {
  return dataset.filter(row =>
    nearlyEqual(row.T_C, T_C) &&
    nearlyEqual(row.ratio_h2o_to_co2, ratio_h2o_to_co2) &&
    nearlyEqual(row.conversion, conversion) &&
    nearlyEqual(row.air_outlet_o2_n2_ratio, air_outlet_o2_n2_ratio)
  );
}

function sortRowsByJ(rows) {
  return [...rows].sort((a, b) => Number(a.j_A_cm2) - Number(b.j_A_cm2));
}

function linearInterp(x, x0, x1, y0, y1) {
  if (x1 === x0) return y0;
  return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}

function interpolateRowAtJ(rows, jQuery) {
  if (!rows || rows.length < 2) {
    throw new Error("Need at least two rows for interpolation.");
  }

  const sorted = sortRowsByJ(rows);
  const jValues = sorted.map(r => Number(r.j_A_cm2));

  let leftRow, rightRow;

  if (jQuery <= jValues[0]) {
    leftRow = sorted[0];
    rightRow = sorted[1];
  } else if (jQuery >= jValues[jValues.length - 1]) {
    leftRow = sorted[sorted.length - 2];
    rightRow = sorted[sorted.length - 1];
  } else {
    for (let i = 0; i < sorted.length - 1; i++) {
      const j0 = Number(sorted[i].j_A_cm2);
      const j1 = Number(sorted[i + 1].j_A_cm2);
      if (j0 <= jQuery && jQuery <= j1) {
        leftRow = sorted[i];
        rightRow = sorted[i + 1];
        break;
      }
    }
  }

  if (!leftRow || !rightRow) {
    throw new Error("Failed to locate interpolation interval.");
  }

  const j0 = Number(leftRow.j_A_cm2);
  const j1 = Number(rightRow.j_A_cm2);

  const out = {};
  const keys = Object.keys(leftRow);

  for (const key of keys) {
    const v0 = leftRow[key];
    const v1 = rightRow[key];
    const n0 = Number(v0);
    const n1 = Number(v1);

    if (Number.isFinite(n0) && Number.isFinite(n1)) {
      out[key] = linearInterp(jQuery, j0, j1, n0, n1);
    } else {
      out[key] = v0;
    }
  }

  out.j_A_cm2 = jQuery;
  return out;
}

function parseIvPoints(text) {
  const lines = text
    .split("\n")
    .map(s => s.trim())
    .filter(Boolean);

  const points = [];

  for (const line of lines) {
    const parts = line.split(",").map(s => s.trim());
    if (parts.length !== 2) {
      throw new Error(`Invalid IV line: "${line}"`);
    }

    const j = Number(parts[0]);
    const v = Number(parts[1]);

    if (!Number.isFinite(j) || !Number.isFinite(v)) {
      throw new Error(`Invalid numeric IV line: "${line}"`);
    }

    points.push({ j_A_cm2: j, V_V: v });
  }

  if (points.length < 2) {
    throw new Error("At least two IV points are required.");
  }

  points.sort((a, b) => a.j_A_cm2 - b.j_A_cm2);

  const unique = [];
  for (const p of points) {
    if (
      unique.length === 0 ||
      Math.abs(p.j_A_cm2 - unique[unique.length - 1].j_A_cm2) > 1e-12
    ) {
      unique.push(p);
    }
  }

  if (unique.length < 2) {
    throw new Error("Need at least two unique j points.");
  }

  return unique;
}

function voltageFromCurrentDensity(ivPoints, jQuery) {
  const pts = [...ivPoints].sort((a, b) => a.j_A_cm2 - b.j_A_cm2);

  if (jQuery <= pts[0].j_A_cm2) {
    return linearInterp(
      jQuery,
      pts[0].j_A_cm2,
      pts[1].j_A_cm2,
      pts[0].V_V,
      pts[1].V_V
    );
  }

  if (jQuery >= pts[pts.length - 1].j_A_cm2) {
    const n = pts.length;
    return linearInterp(
      jQuery,
      pts[n - 2].j_A_cm2,
      pts[n - 1].j_A_cm2,
      pts[n - 2].V_V,
      pts[n - 1].V_V
    );
  }

  for (let i = 0; i < pts.length - 1; i++) {
    const p0 = pts[i];
    const p1 = pts[i + 1];
    if (p0.j_A_cm2 <= jQuery && jQuery <= p1.j_A_cm2) {
      return linearInterp(
        jQuery,
        p0.j_A_cm2,
        p1.j_A_cm2,
        p0.V_V,
        p1.V_V
      );
    }
  }

  throw new Error("Failed to interpolate IV voltage.");
}

function stackPowerW(j_A_cm2, V_V, n_cells, area_cm2) {
  return n_cells * area_cm2 * j_A_cm2 * V_V;
}

function computeResidualAtJ(rows, ivPoints, j_A_cm2, n_cells, area_cm2) {
  const thermoRow = interpolateRowAtJ(rows, j_A_cm2);
  const V_V = voltageFromCurrentDensity(ivPoints, j_A_cm2);
  const P_stack_W = stackPowerW(j_A_cm2, V_V, n_cells, area_cm2);
  const residual_W = Number(thermoRow.thermo_term_W) - P_stack_W;

  return {
    j_A_cm2,
    V_V,
    P_stack_W,
    thermo_term_W: Number(thermoRow.thermo_term_W),
    residual_W,
    thermoRow,
  };
}

function findBracketForRoot(rows, ivPoints, n_cells, area_cm2) {
  const sorted = sortRowsByJ(rows);
  const jMin = Number(sorted[0].j_A_cm2);
  const jMax = Number(sorted[sorted.length - 1].j_A_cm2);

  const nScan = 200;
  let prev = null;

  for (let i = 0; i <= nScan; i++) {
    const j = jMin + (jMax - jMin) * i / nScan;
    const cur = computeResidualAtJ(rows, ivPoints, j, n_cells, area_cm2);

    if (prev && prev.residual_W === 0) {
      return [prev.j_A_cm2, prev.j_A_cm2];
    }

    if (prev && prev.residual_W * cur.residual_W < 0) {
      return [prev.j_A_cm2, cur.j_A_cm2];
    }

    prev = cur;
  }

  return null;
}

function bisectRoot(rows, ivPoints, n_cells, area_cm2, jLeft, jRight, tol = 1e-6, maxIter = 100) {
  let left = jLeft;
  let right = jRight;

  let fLeft = computeResidualAtJ(rows, ivPoints, left, n_cells, area_cm2).residual_W;
  let fRight = computeResidualAtJ(rows, ivPoints, right, n_cells, area_cm2).residual_W;

  if (Math.abs(fLeft) < tol) {
    return computeResidualAtJ(rows, ivPoints, left, n_cells, area_cm2);
  }
  if (Math.abs(fRight) < tol) {
    return computeResidualAtJ(rows, ivPoints, right, n_cells, area_cm2);
  }
  if (fLeft * fRight > 0) {
    throw new Error("Root is not bracketed.");
  }

  for (let iter = 0; iter < maxIter; iter++) {
    const mid = 0.5 * (left + right);
    const fMidObj = computeResidualAtJ(rows, ivPoints, mid, n_cells, area_cm2);
    const fMid = fMidObj.residual_W;

    if (Math.abs(fMid) < tol || Math.abs(right - left) < tol) {
      return fMidObj;
    }

    if (fLeft * fMid < 0) {
      right = mid;
      fRight = fMid;
    } else {
      left = mid;
      fLeft = fMid;
    }
  }

  return computeResidualAtJ(rows, ivPoints, 0.5 * (left + right), n_cells, area_cm2);
}

function findTNPoint(rows, ivPoints, n_cells, area_cm2) {
  const bracket = findBracketForRoot(rows, ivPoints, n_cells, area_cm2);
  if (!bracket) {
    throw new Error("No TN root found within dataset j range.");
  }

  const [jLeft, jRight] = bracket;
  return bisectRoot(rows, ivPoints, n_cells, area_cm2, jLeft, jRight);
}

function getInputValue(id) {
  return Number(document.getElementById(id).value);
}

function getInputText(id) {
  return document.getElementById(id).value;
}

function render(obj) {
  document.getElementById("output").textContent = JSON.stringify(obj, null, 2);
}

let DATASET = [];

function getSelectedCondition() {
  return {
    T_C: getInputValue("T_C"),
    ratio_h2o_to_co2: getInputValue("ratio"),
    conversion: getInputValue("conversion"),
    air_outlet_o2_n2_ratio: getInputValue("air_ratio"),
  };
}

function getStackSetting() {
  return {
    n_cells: getInputValue("n_cells"),
    area_cm2: getInputValue("area_cm2"),
  };
}

async function init() {
  try {
    DATASET = await loadThermoDataset("./thermo_dataset.json");
    render({ status: "dataset loaded", n_records: DATASET.length });
  } catch (err) {
    render({ error: String(err) });
    return;
  }

  document.getElementById("runBtn").addEventListener("click", () => {
    try {
      const condition = getSelectedCondition();
      const jQuery = getInputValue("j_query");

      const selectedRows = filterRowsByCondition(DATASET, condition);
      if (selectedRows.length === 0) {
        throw new Error("No matching rows found for the selected condition.");
      }

      const interpolated = interpolateRowAtJ(selectedRows, jQuery);

      render({
        mode: "interpolation only",
        matched_rows: selectedRows.length,
        interpolated_row: interpolated,
      });
    } catch (err) {
      render({ error: String(err) });
    }
  });

  document.getElementById("tnBtn").addEventListener("click", () => {
    try {
      const condition = getSelectedCondition();
      const { n_cells, area_cm2 } = getStackSetting();
      const ivPoints = parseIvPoints(getInputText("iv_points"));

      const selectedRows = filterRowsByCondition(DATASET, condition);
      if (selectedRows.length === 0) {
        throw new Error("No matching rows found for the selected condition.");
      }

      const tn = findTNPoint(selectedRows, ivPoints, n_cells, area_cm2);

      const H2_mol_s = Number(tn.thermoRow.fuel_out_n_H2_mol_s || 0);
      const CO_mol_s = Number(tn.thermoRow.fuel_out_n_CO_mol_s || 0);

      render({
        mode: "TN solver",
        matched_rows: selectedRows.length,
        TN_j_A_cm2: tn.j_A_cm2,
        TN_V_V: tn.V_V,
        P_stack_W: tn.P_stack_W,
        thermo_term_W: tn.thermo_term_W,
        residual_W: tn.residual_W,
        fuel_out_n_H2_mol_s: H2_mol_s,
        fuel_out_n_CO_mol_s: CO_mol_s,
        thermo_row_at_TN: tn.thermoRow,
      });
    } catch (err) {
      render({ error: String(err) });
    }
  });
}

init();