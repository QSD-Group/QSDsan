(function () {
  function getColumnIndex(table, name) {
    const headers = Array.from(table.querySelectorAll("thead th"));
    return headers.findIndex((header) => header.textContent.trim() === name);
  }

  function initUnitOperationFilters() {
    const table = document.querySelector("table.unit-operation-table");
    const categoryFilter = document.getElementById("unit-operation-category-filter");
    const abstractFilter = document.getElementById("unit-operation-abstract-filter");
    const count = document.getElementById("unit-operation-filter-count");

    if (!table || !categoryFilter || !abstractFilter) return;

    const categoryIndex = getColumnIndex(table, "Category");
    const abstractIndex = getColumnIndex(table, "Abstract");
    if (categoryIndex < 0 || abstractIndex < 0) return;

    const rows = Array.from(table.querySelectorAll("tbody tr"));

    function applyFilters() {
      const category = categoryFilter.value;
      const abstractValue = abstractFilter.value;
      let visible = 0;

      rows.forEach((row) => {
        const cells = row.querySelectorAll("td");
        const rowCategory = cells[categoryIndex]?.textContent.trim() || "";
        const rowAbstract = cells[abstractIndex]?.textContent.trim() || "";
        const show =
          (!category || rowCategory === category) &&
          (!abstractValue || rowAbstract === abstractValue);

        row.hidden = !show;
        if (show) visible += 1;
      });

      if (count) count.textContent = visible + " shown";
    }

    categoryFilter.addEventListener("change", applyFilters);
    abstractFilter.addEventListener("change", applyFilters);
    applyFilters();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initUnitOperationFilters);
  } else {
    initUnitOperationFilters();
  }
})();
