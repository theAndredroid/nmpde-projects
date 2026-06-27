#!/usr/bin/env pvpython
"""
Script to generate clip plane visualizations of activation time for model sweeps.
Applies a ParaView Plane Clip centered at (10, 3.5, 1.5) with normal vector defined by
the cross product of (20, 3.5, 1.5) and (-20, 3.5, 0), colored using the 'Jet' colormap.
Camera is positioned starting from corner (0, 0, 3) looking towards the center (10, 3.5, 1.5).
All model outputs share a unified color scale spanning the global min and max activation times.
"""

import sys
import os
os.environ['DISPLAY'] = ''  # Force headless offscreen rendering to prevent window flashing
import glob
import re
import argparse
import numpy as np

try:
    from paraview.simple import *
    import paraview.servermanager as sm
    HAS_PARAVIEW = True
except ImportError:
    HAS_PARAVIEW = False

def find_latest_models_dir(base_dir="build"):
    candidates = glob.glob(os.path.join(base_dir, "*_models"))
    if not candidates:
        candidates = glob.glob("*_models")
    if not candidates:
        return None
    return max(candidates, key=os.path.getmtime)

def get_file_range(pvtu_path):
    """Extract scalar range for Activation Time from a VTU/PVTU file."""
    if pvtu_path.endswith('.pvtu'):
        reader = XMLPartitionedUnstructuredGridReader(FileName=[pvtu_path])
    else:
        reader = XMLUnstructuredGridReader(FileName=[pvtu_path])
        
    reader.UpdatePipeline()
    fetched = sm.Fetch(reader)
    arr = fetched.GetPointData().GetArray('Activation Time')
    rng = arr.GetRange() if arr else None
    Delete(reader)
    return rng

def render_model_clip(pvtu_path, output_png, min_val=None, max_val=None):
    """Render activation time on a clipped plane with Jet colormap and fixed color scale."""
    if not HAS_PARAVIEW:
        raise RuntimeError("ParaView python modules are not available.")

    # Load reader
    if pvtu_path.endswith('.pvtu'):
        reader = XMLPartitionedUnstructuredGridReader(FileName=[pvtu_path])
    else:
        reader = XMLUnstructuredGridReader(FileName=[pvtu_path])

    # Define Clip plane parameters
    # Origin: (10, 3.5, 1.5)
    origin = [10.0, 3.5, 1.5]
    
    # Normal: Cross product of (20, 3.5, 1.5) and (-20, 3.5, 0)
    v1 = np.array([20.0, 7, 3])
    v2 = np.array([-20.0, 7, 0.0])
    normal = list(np.cross(v1, v2))

    clip = Clip(Input=reader)
    clip.ClipType = 'Plane'
    clip.ClipType.Origin = origin
    clip.ClipType.Normal = normal

    # Create render view
    renderView = CreateView('RenderView')
    renderView.ViewSize = [1600, 1000]
    renderView.Background = [0.15, 0.15, 0.15]  # Dark background for contrast

    display = Show(clip, renderView)
    display.Representation = 'Surface'
    ColorBy(display, ('POINTS', 'Activation Time'))

    # Apply Jet color preset and unify color scale across models
    lut = GetColorTransferFunction('ActivationTime')
    lut.ApplyPreset('Jet', True)
    if min_val is not None and max_val is not None:
        lut.RescaleTransferFunction(min_val, max_val)
    display.LookupTable = lut

    # Show colorbar scalar bar
    scalarBar = GetScalarBar(lut, renderView)
    scalarBar.Title = 'Activation Time [ms]'
    scalarBar.ComponentTitle = ''
    scalarBar.TitleColor = [0.0, 0.0, 0.0]
    scalarBar.LabelColor = [0.0, 0.0, 0.0]
    scalarBar.Visibility = 1

    # Position camera at corner (0, 0, 3) looking towards domain center (10, 3.5, 1.5)
    renderView.CameraPosition = [-17.0525, -16.2422, 25.7423]
    renderView.CameraFocalPoint = origin
    renderView.CameraViewUp = [0.444443, 0.384503, 0.809091]

    Render()
    SaveScreenshot(output_png, renderView, ImageResolution=[1600, 1000], TransparentBackground=1)
    
    # Cleanup ParaView objects
    Delete(renderView)
    Delete(clip)
    Delete(reader)

def main():
    parser = argparse.ArgumentParser(description="Generate clipped activation time PNG plots for model sweeps.")
    parser.add_argument("models_dir", nargs="?", default=None, help="Path to the models sweep directory")
    args = parser.parse_args()

    models_dir = args.models_dir
    if not models_dir:
        models_dir = find_latest_models_dir()

    if not models_dir or not os.path.exists(models_dir):
        print(f"Error: Could not locate models sweep directory: {models_dir}")
        sys.exit(1)

    print(f"Processing model sweep directory: {models_dir}")

    subdirs = sorted(glob.glob(os.path.join(models_dir, "*")))
    tasks = []

    # 1. Collect master files only (.pvtu files preferred, skipping sub-partition pieces like _0.00.vtu)
    for sd in subdirs:
        if not os.path.isdir(sd):
            continue
        pvtu_files = glob.glob(os.path.join(sd, "*.pvtu"))
        if not pvtu_files:
            all_vtu = glob.glob(os.path.join(sd, "*.vtu"))
            pvtu_files = [f for f in all_vtu if not re.search(r'\.[0-9]{2,}\.vtu$', f)]
            
        for pvtu in sorted(list(set(pvtu_files))):
            tasks.append((sd, pvtu))

    if not tasks:
        print("Warning: No VTU/PVTU master files found to process.")
        sys.exit(0)

    # 2. Compute Global Min and Max across all model datasets
    print("Computing global scalar range for unified colormap scale...")
    all_mins, all_maxs = [], []
    for sd, pvtu in tasks:
        rng = get_file_range(pvtu)
        if rng:
            all_mins.append(rng[0])
            all_maxs.append(rng[1])

    global_min = min(all_mins) if all_mins else None
    global_max = max(all_maxs) if all_maxs else None
    if global_min is not None and global_max is not None:
        print(f"  [+] Unified Color Scale Range: [{global_min:.4f}, {global_max:.4f}] ms")

    # 3. Render master screenshots with unified color scale
    for sd, pvtu in tasks:
        base_name = os.path.splitext(os.path.basename(pvtu))[0]
        output_png = os.path.join(models_dir, f"{base_name}_clip.png")
        print(f"  [+] Rendering: {pvtu}")
        try:
            render_model_clip(pvtu, output_png, min_val=global_min, max_val=global_max)
            print(f"      -> Saved PNG: {output_png}")
        except Exception as e:
            print(f"      -> Error rendering {pvtu}: {e}")

if __name__ == "__main__":
    main()
