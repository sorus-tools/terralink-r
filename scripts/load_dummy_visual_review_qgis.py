import csv
import json
import os
import sys

from qgis.PyQt.QtGui import QColor
from qgis.core import (
    QgsFillSymbol,
    QgsLayerTreeGroup,
    QgsLayerTreeLayer,
    QgsPalettedRasterRenderer,
    QgsProject,
    QgsRasterLayer,
    QgsSingleSymbolRenderer,
    QgsVectorLayer,
)


PLUGIN_PARENT = "/Users/benbishop/Library/Application Support/QGIS/QGIS3/profiles/default/python/plugins"
if PLUGIN_PARENT not in sys.path:
    sys.path.insert(0, PLUGIN_PARENT)

from TerraLink_v1_7.analysis_raster import run_raster_analysis
from TerraLink_v1_7.analysis_vector import run_vector_analysis


STRATEGIES = [
    "most_connected_habitat",
    "largest_single_network",
    "reachable_habitat_advanced",
    "landscape_fluidity",
]


def _as_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _remove_group(root: QgsLayerTreeGroup, name: str) -> None:
    for child in list(root.children()):
        if isinstance(child, QgsLayerTreeGroup) and child.name() == name:
            root.removeChildNode(child)


def _ensure_group(parent: QgsLayerTreeGroup, name: str) -> QgsLayerTreeGroup:
    for child in parent.children():
        if isinstance(child, QgsLayerTreeGroup) and child.name() == name:
            return child
    group = parent.addGroup(name)
    try:
        group.setExpanded(False)
    except Exception:
        pass
    return group


def _load_vector_layer(uri: str, name: str, group: QgsLayerTreeGroup, fill: str, outline: str, visible: bool = False) -> QgsVectorLayer:
    layer = QgsVectorLayer(uri, name, "ogr")
    if not layer.isValid():
      raise RuntimeError(f"Failed to load vector layer: {uri}")
    symbol = QgsFillSymbol.createSimple(
        {
            "color": fill,
            "outline_color": outline,
            "outline_width": "0.4",
        }
    )
    layer.setRenderer(QgsSingleSymbolRenderer(symbol))
    QgsProject.instance().addMapLayer(layer, False)
    node = group.addLayer(layer)
    node.setItemVisibilityChecked(bool(visible))
    return layer


def _load_raster_layer(path: str, name: str, group: QgsLayerTreeGroup, visible: bool = False) -> QgsRasterLayer:
    layer = QgsRasterLayer(path, name)
    if not layer.isValid():
        raise RuntimeError(f"Failed to load raster layer: {path}")
    provider = layer.dataProvider()
    classes = [
        QgsPalettedRasterRenderer.Class(0, QColor(255, 255, 255, 0), "background"),
        QgsPalettedRasterRenderer.Class(1, QColor("#3A7D44"), "habitat"),
        QgsPalettedRasterRenderer.Class(9, QColor("#1F1F1F"), "impassable"),
    ]
    try:
        layer.setRenderer(QgsPalettedRasterRenderer(provider, 1, classes))
    except Exception:
        pass
    QgsProject.instance().addMapLayer(layer, False)
    node = group.addLayer(layer)
    node.setItemVisibilityChecked(bool(visible))
    return layer


def _run_group_label(strategy: str, stats: dict) -> str:
    corridors = int((stats or {}).get("corridors_used", 0) or 0)
    budget = float((stats or {}).get("budget_used", 0.0) or 0.0)
    return f"{strategy} [c={corridors}, b={budget:.2f}]"


def _record_stats(stats: dict) -> dict:
    return {
        "corridors_used": int((stats or {}).get("corridors_used", 0) or 0),
        "budget_used": float(
            (stats or {}).get("budget_used_ha", (stats or {}).get("budget_used_display", (stats or {}).get("budget_used", 0.0)))
            or 0.0
        ),
        "total_connected": float(
            (stats or {}).get("total_connected_area_ha", (stats or {}).get("total_connected_area_display", 0.0)) or 0.0
        ),
        "largest_network": float(
            (stats or {}).get("largest_group_area_ha", (stats or {}).get("largest_group_area_display", 0.0)) or 0.0
        ),
        "habitat_availability_post": float((stats or {}).get("habitat_availability_after", 0.0) or 0.0),
        "mean_reachable_area_post": float((stats or {}).get("mean_reachable_area", 0.0) or 0.0),
        "largest_reachable_habitat_cluster_post": float(
            (stats or {}).get("largest_reachable_habitat_cluster", (stats or {}).get("largest_reachable_habitat_cluster_display", 0.0))
            or 0.0
        ),
        "mean_effective_resistance_post": float((stats or {}).get("mean_effective_resistance_post_exact", 0.0) or 0.0),
        "landscape_fluidity_post": float((stats or {}).get("landscape_fluidity_post_exact", 0.0) or 0.0),
    }


def run_visual_review(manifest_path: str) -> dict:
    with open(manifest_path, "r", encoding="utf-8") as fh:
        manifest = json.load(fh)

    root = QgsProject.instance().layerTreeRoot()
    group_name = manifest["group_name"]
    _remove_group(root, group_name)
    top_group = root.addGroup(group_name)

    inputs_group = _ensure_group(top_group, "Inputs")
    vector_input = _load_vector_layer(
        f"{manifest['inputs']['vector_patches_path']}|layername={manifest['inputs']['vector_patches_layer']}",
        "dummy vector patches",
        inputs_group,
        "110,168,84,90",
        "#295135",
        visible=True,
    )
    obstacle_input = _load_vector_layer(
        f"{manifest['inputs']['vector_obstacles_path']}|layername={manifest['inputs']['vector_obstacles_layer']}",
        "dummy vector impassable",
        inputs_group,
        "40,40,40,160",
        "#111111",
        visible=True,
    )
    raster_input = _load_raster_layer(
        manifest["inputs"]["raster_path"],
        "dummy raster",
        inputs_group,
        visible=False,
    )

    r_engine_group = _ensure_group(top_group, "R Results")
    qgis_engine_group = _ensure_group(top_group, "QGIS Results")

    for run in manifest.get("r_runs", []):
        parent = _ensure_group(r_engine_group, run["input_type"].title())
        parent = _ensure_group(parent, run["obstacle_mode"].replace("_", " ").title())
        run_group = _ensure_group(parent, _run_group_label(run["strategy"], run["stats"]))
        _load_vector_layer(
            f"{run['gpkg_path']}|layername={run['corridor_layer']}",
            "corridors",
            run_group,
            "67,113,201,90",
            "#254C8A",
            visible=False,
        )
        _load_vector_layer(
            f"{run['gpkg_path']}|layername={run['network_layer']}",
            "networks",
            run_group,
            "67,113,201,45",
            "#16325C",
            visible=False,
        )
        if run["input_type"] == "raster":
            if run.get("corridor_raster_path"):
                _load_raster_layer(run["corridor_raster_path"], "corridor raster", run_group, visible=False)
            if run.get("contiguous_raster_path"):
                _load_raster_layer(run["contiguous_raster_path"], "contiguous raster", run_group, visible=False)

    qgis_records = []
    comparison_rows = []

    vector_params_base = {
        "min_corridor_width": manifest["params"]["vector"]["min_corridor_width"],
        "min_patch_size": manifest["params"]["vector"]["min_patch_size"],
        "budget_area": manifest["params"]["vector"]["budget"],
        "max_search_distance": manifest["params"]["vector"]["max_search_distance"],
        "unit_system": manifest["params"]["vector"]["units"],
        "grid_resolution": 30.0,
        "species_dispersal_distance_analysis": manifest["params"]["vector"]["species_dispersal_distance"],
        "species_dispersal_kernel": manifest["params"]["vector"]["species_dispersal_kernel"],
        "min_patch_area_for_species_analysis": 0.0,
        "patch_quality_weight_field": manifest["params"]["vector"]["patch_quality_field"],
        "patch_area_scaling": manifest["params"]["vector"]["patch_area_scaling"],
        "add_to_project": False,
    }
    raster_params_base = {
        "patch_connectivity": 8,
        "patch_mode": "value",
        "patch_values": _as_list(manifest["params"]["raster"]["patch_values"]),
        "obstacle_mode": "value",
        "obstacle_values": _as_list(manifest["params"]["raster"]["obstacle_values"]),
        "min_patch_size": manifest["params"]["raster"]["min_patch_size"],
        "budget_pixels": manifest["params"]["raster"]["budget"],
        "max_search_distance": manifest["params"]["raster"]["max_search_distance"],
        "min_corridor_width": manifest["params"]["raster"]["min_corridor_width"],
        "allow_bottlenecks": False,
        "species_dispersal_distance_analysis": manifest["params"]["raster"]["species_dispersal_distance"],
        "species_dispersal_kernel": manifest["params"]["raster"]["species_dispersal_kernel"],
        "min_patch_area_for_species_analysis": 0.0,
        "patch_area_scaling": manifest["params"]["raster"]["patch_area_scaling"],
        "add_to_project": False,
    }

    for input_type in ("vector", "raster"):
        for obstacle_mode in ("no_obstacles", "with_obstacles"):
            parent = _ensure_group(qgis_engine_group, input_type.title())
            parent = _ensure_group(parent, obstacle_mode.replace("_", " ").title())
            for strategy in STRATEGIES:
                run_key = f"qgis__{input_type}__{obstacle_mode}__{strategy}"
                output_dir = os.path.join(manifest["root_dir"], "qgis_runs", run_key)
                os.makedirs(output_dir, exist_ok=True)

                if input_type == "vector":
                    raw_params = dict(vector_params_base)
                    raw_params["output_name"] = f"{run_key}.gpkg"
                    raw_params["obstacle_enabled"] = obstacle_mode == "with_obstacles"
                    raw_params["obstacle_layer_ids"] = [obstacle_input.id()] if obstacle_mode == "with_obstacles" else []
                    result = run_vector_analysis(
                        vector_input,
                        output_dir,
                        raw_params,
                        strategy=strategy,
                        temporary=False,
                        iface=None,
                    )
                else:
                    raw_params = dict(raster_params_base)
                    raw_params["obstacle_enabled"] = obstacle_mode == "with_obstacles"
                    if obstacle_mode != "with_obstacles":
                        raw_params["obstacle_values"] = []
                    result = run_raster_analysis(
                        raster_input,
                        output_dir,
                        raw_params,
                        strategy=strategy,
                        temporary=False,
                        iface=None,
                    )

                item = result[0]
                stats = item.get("stats") or {}
                gpkg_path = item["output_path"]
                corridors_layer_name = item["layer_name"]
                networks_layer_name = "Contiguous Areas"
                stats_out = _record_stats(stats)
                run_group = _ensure_group(parent, _run_group_label(strategy, stats_out))
                _load_vector_layer(
                    f"{gpkg_path}|layername={corridors_layer_name}",
                    "corridors",
                    run_group,
                    "214,108,44,90",
                    "#8A4317",
                    visible=False,
                )
                _load_vector_layer(
                    f"{gpkg_path}|layername={networks_layer_name}",
                    "networks",
                    run_group,
                    "214,108,44,45",
                    "#6B310F",
                    visible=False,
                )

                qgis_record = {
                    "engine": "QGIS",
                    "input_type": input_type,
                    "obstacle_mode": obstacle_mode,
                    "strategy": strategy,
                    "run_key": run_key,
                    "display_label": f"QGIS | {input_type} | {obstacle_mode} | {strategy}",
                    "gpkg_path": gpkg_path,
                    "corridor_layer": corridors_layer_name,
                    "network_layer": networks_layer_name,
                    "stats": stats_out,
                }
                qgis_records.append(qgis_record)

    qgis_summary_path = os.path.join(manifest["root_dir"], "qgis_run_summary.json")
    with open(qgis_summary_path, "w", encoding="utf-8") as fh:
        json.dump(qgis_records, fh, indent=2, sort_keys=True)

    qgis_index = {
        (row["input_type"], row["obstacle_mode"], row["strategy"]): row
        for row in qgis_records
    }
    for row in manifest.get("r_runs", []):
        key = (row["input_type"], row["obstacle_mode"], row["strategy"])
        qrow = qgis_index.get(key)
        if qrow is None:
            continue
        comparison_rows.append(
            {
                "input_type": row["input_type"],
                "obstacle_mode": row["obstacle_mode"],
                "strategy": row["strategy"],
                "r_corridors_used": row["stats"]["corridors_used"],
                "qgis_corridors_used": qrow["stats"]["corridors_used"],
                "r_budget_used": row["stats"]["budget_used"],
                "qgis_budget_used": qrow["stats"]["budget_used"],
                "r_total_connected": row["stats"]["total_connected"],
                "qgis_total_connected": qrow["stats"]["total_connected"],
                "r_largest_network": row["stats"]["largest_network"],
                "qgis_largest_network": qrow["stats"]["largest_network"],
                "r_landscape_fluidity_post": row["stats"]["landscape_fluidity_post"],
                "qgis_landscape_fluidity_post": qrow["stats"]["landscape_fluidity_post"],
            }
        )

    comparison_csv_path = os.path.join(manifest["root_dir"], "comparison_summary.csv")
    with open(comparison_csv_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "input_type",
                "obstacle_mode",
                "strategy",
                "r_corridors_used",
                "qgis_corridors_used",
                "r_budget_used",
                "qgis_budget_used",
                "r_total_connected",
                "qgis_total_connected",
                "r_largest_network",
                "qgis_largest_network",
                "r_landscape_fluidity_post",
                "qgis_landscape_fluidity_post",
            ],
        )
        writer.writeheader()
        writer.writerows(comparison_rows)

    try:
        if "iface" in globals() and iface is not None:
            iface.mapCanvas().setExtent(vector_input.extent())
            iface.mapCanvas().refresh()
    except Exception:
        pass

    return {
        "manifest_path": manifest_path,
        "qgis_summary_path": qgis_summary_path,
        "comparison_csv_path": comparison_csv_path,
        "group_name": group_name,
        "qgis_runs": len(qgis_records),
    }


if __name__ == "__main__":
    manifest_path = globals().get("MANIFEST_PATH")
    if not manifest_path:
        raise RuntimeError("Set MANIFEST_PATH before executing this script.")
    result = run_visual_review(manifest_path)
    print(json.dumps(result, sort_keys=True))
