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

let DATASET = [];

function getInputValue(id) {
  return Number(document.getElementById(id).value);
}

function render(obj) {
  document.getElementById("output").textContent = JSON.stringify(obj, null, 2);
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
      const condition = {
        T_C: getInputValue("T_C"),
        ratio_h2o_to_co2: getInputValue("ratio"),
        conversion: getInputValue("conversion"),
        air_outlet_o2_n2_ratio: getInputValue("air_ratio"),
      };

      const jQuery = getInputValue("j_query");

      const selectedRows = filterRowsByCondition(DATASET, condition);

      if (selectedRows.length === 0) {
        throw new Error("No matching rows found for the selected condition.");
      }

      const interpolated = interpolateRowAtJ(selectedRows, jQuery);

      render({
        matched_rows: selectedRows.length,
        interpolated_row: interpolated,
      });
    } catch (err) {
      render({ error: String(err) });
    }
  });
}

init();