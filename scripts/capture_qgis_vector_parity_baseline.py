#!/usr/bin/env python3

import argparse
import json
import subprocess
import textwrap
from pathlib import Path


QGIS_SEND = "/Users/benbishop/.codex/skills/QGIS-connect/scripts/qgis_send.py"
PLUGIN_PARENT = "/Users/benbishop/Library/Application Support/QGIS/QGIS3/profiles/default/python/plugins"
VECTOR_FIXTURE = "/Users/benbishop/projects/terralink/inst/extdata/synthetic_patches.gpkg|layername=patches"
OBSTACLE_FIXTURE = "/Users/benbishop/projects/terralink/inst/extdata/synthetic_impassable.gpkg|layername=impassable"


def build_qgis_code() -> str:
    return textwrap.dedent(
        f"""
        import json, sys, tempfile
        from qgis.core import QgsVectorLayer, QgsProject
        parent_dir = r"{PLUGIN_PARENT}"
        if parent_dir not in sys.path:
            sys.path.insert(0, parent_dir)
        from TerraLink_v1_7.analysis_vector import run_vector_analysis

        vec_path = r"{VECTOR_FIXTURE}"
        obs_path = r"{OBSTACLE_FIXTURE}"
        strategies = ["most_connected_habitat", "largest_single_network", "reachable_habitat_advanced", "landscape_fluidity"]
        vec = QgsVectorLayer(vec_path, "synthetic_vector_fixture", "ogr")
        obs = QgsVectorLayer(obs_path, "synthetic_obstacle_fixture", "ogr")
        project = QgsProject.instance()
        project.addMapLayer(obs, False)
        params = {{
            "min_corridor_width": 60.0,
            "min_patch_size": 0.01,
            "budget_area": 18.0,
            "max_search_distance": 900.0,
            "unit_system": "metric",
            "obstacle_enabled": True,
            "obstacle_layer_ids": [obs.id()],
            "obstacle_layer_id": obs.id(),
            "grid_resolution": 30.0,
            "species_dispersal_distance_analysis": 800.0,
            "species_dispersal_kernel": "exponential",
            "min_patch_area_for_species_analysis": 0.0,
            "patch_quality_weight_field": "quality",
            "patch_area_scaling": "sqrt",
            "add_to_project": False,
        }}
        out = {{}}
        try:
            for strategy in strategies:
                tmpdir = tempfile.mkdtemp(prefix=f"qgis_vector_{{strategy}}_")
                result = run_vector_analysis(vec, tmpdir, dict(params), strategy=strategy, temporary=False, iface=None)
                item = result[0]
                stats = item.get("stats") or {{}}
                layer_name = item.get("layer_name")
                output_path = item.get("output_path")
                corridor_layer = QgsVectorLayer(f"{{output_path}}|layername={{layer_name}}", f"corr_{{strategy}}", "ogr")
                rows = []
                for feat in corridor_layer.getFeatures():
                    rows.append({{
                        "patch1": int(feat["patch1"]),
                        "patch2": int(feat["patch2"]),
                        "corridor_area_ha": float(feat["corridor_area_ha"]),
                        "patches_area_ha": float(feat["patches_area_ha"]),
                        "network_area_ha": float(feat["network_area_ha"]),
                        "efficiency": float(feat["efficiency"]),
                        "redundant": str(feat["redundant"]),
                    }})
                rows.sort(key=lambda r: (r["patch1"], r["patch2"]))
                out[strategy] = {{
                    "stats": {{
                        k: stats.get(k) for k in [
                            "corridors_used",
                            "budget_used_ha",
                            "total_connected_area_ha",
                            "largest_group_area_ha",
                            "habitat_availability_after",
                            "mean_reachable_area",
                            "largest_reachable_habitat_cluster",
                            "mean_effective_resistance_pre_exact",
                            "mean_effective_resistance_post_exact",
                            "landscape_fluidity_pre_exact",
                            "landscape_fluidity_post_exact",
                        ]
                    }},
                    "corridors": rows,
                }}
            print(json.dumps(out, sort_keys=True, indent=2))
        finally:
            project.removeMapLayer(obs.id())
        """
    ).strip()


def main() -> int:
    parser = argparse.ArgumentParser(description="Capture the current QGIS vector parity baseline for TerraLink.")
    parser.add_argument("--output", type=Path, help="Optional output JSON path.")
    args = parser.parse_args()

    proc = subprocess.run(
        ["python3", QGIS_SEND, "--code", build_qgis_code()],
        check=True,
        capture_output=True,
        text=True,
    )
    response = json.loads(proc.stdout)
    if not response.get("ok"):
        raise SystemExit(response)
    payload = response.get("stdout", "").strip()
    if args.output:
        args.output.write_text(payload + "\n", encoding="utf-8")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
