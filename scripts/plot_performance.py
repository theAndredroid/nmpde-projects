#!/usr/bin/env python3
import sys
import os
import re

def parse_results(results_file):
    data = []
    pattern = re.compile(r'^(\S+)\s+-\s+(.*?):\s+([0-9.]+)\s+seconds')
    with open(results_file, 'r') as f:
        for line in f:
            match = pattern.match(line.strip())
            if match:
                job_id = match.group(1)
                config = match.group(2).strip()
                # Clean up multiple spaces
                config = re.sub(r'\s+', ' ', config)
                time_s = float(match.group(3))
                data.append((job_id, config, time_s))
    return data

def parse_iterations(out_file):
    iterations = []
    pattern = re.compile(r'Timestep\s+\d+.*:\s*(\d+)\s+\w+\s+iterations')
    if not os.path.isfile(out_file):
        return 0.0
    with open(out_file, 'r') as f:
        for line in f:
            match = pattern.search(line)
            if match:
                iterations.append(int(match.group(1)))
    return sum(iterations) / len(iterations) if iterations else 0.0

def text_plot(data):
    if not data:
        print("No data to display.")
        return
    print("\n=== Text-Based Performance Plot ===")
    max_len = max(len(d[1]) for d in data)
    
    # Times scale
    max_time = max(d[2] for d in data)
    scale_time = 25.0 / max_time if max_time > 0 else 1.0
    
    # Iters scale
    max_iter = max(d[3] for d in data)
    scale_iter = 25.0 / max_iter if max_iter > 0 else 1.0
    
    for job_id, config, time_s, avg_iter in data:
        bar_time = '#' * int(time_s * scale_time)
        bar_iter = '*' * int(avg_iter * scale_iter)
        print(f"{config:<{max_len}} | Time: {bar_time:<25} ({time_s:.2f}m) | Iters: {bar_iter:<25} ({avg_iter:.1f})")
    print("===================================\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: plot_performance.py <results_directory>")
        sys.exit(1)
        
    directory = sys.argv[1]
    results_file = os.path.join(directory, "solver_performance.txt")
    if not os.path.isfile(results_file):
        print(f"Error: Results file '{results_file}' not found.")
        sys.exit(1)
        
    raw_data = parse_results(results_file)
    if not raw_data:
        print(f"Error: No valid result data found in '{results_file}'.")
        sys.exit(1)
        
    # Process iterations
    data = []
    for job_id, config, time_s in raw_data:
        out_file = os.path.join(directory, f"{job_id}.out")
        avg_iter = parse_iterations(out_file)
        data.append((job_id, config, time_s / 60.0, avg_iter))
        
    # Attempt to plot with matplotlib
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        configs = [d[1] for d in data]
        times = [d[2] for d in data]
        iters = [d[3] for d in data]
        
        # Modern styling
        plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
        
        # Slate/teal palette
        color_palette = ['#2b5c8f', '#4682b4', '#5f9ea0', '#66c2a5', '#3288bd', '#5e4fa2']
        colors = [color_palette[i % len(color_palette)] for i in range(len(data))]
        
        # Plot 1: Execution Time (Bar Chart)
        fig1, ax1 = plt.subplots(figsize=(10, max(5, len(data) * 0.8)))
        bars1 = ax1.barh(configs, times, color=colors, edgecolor='none', height=0.6)
        max_time = max(times) if times else 1.0
        for bar in bars1:
            width = bar.get_width()
            ax1.text(width + (max_time * 0.015), bar.get_y() + bar.get_height()/2,
                     f'{width:.2f}m',
                     va='center', ha='left', fontsize=10, fontweight='bold', color='#2c3e50')
            
        ax1.set_xlabel('Execution Time (minutes)', fontsize=12, fontweight='bold', labelpad=10)
        ax1.set_ylabel('Configuration', fontsize=12, fontweight='bold', labelpad=10)
        ax1.set_title('Solver & Preconditioner Sweep Execution Times', fontsize=14, fontweight='bold', pad=15)
        ax1.grid(True, linestyle='--', alpha=0.5, axis='x')
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_color('#bdc3c7')
        ax1.spines['bottom'].set_color('#bdc3c7')
        plt.tight_layout()
        output_image1 = os.path.join(directory, "execution_times.png")
        plt.savefig(output_image1, dpi=300)
        plt.close(fig1)
        print(f"Successfully generated execution times plot: {output_image1}")
        
        # Plot 2: Average Iterations (Bar Chart)
        fig2, ax2 = plt.subplots(figsize=(10, max(5, len(data) * 0.8)))
        bars2 = ax2.barh(configs, iters, color=colors, edgecolor='none', height=0.6)
        max_iter = max(iters) if iters else 1.0
        for bar in bars2:
            width = bar.get_width()
            ax2.text(width + (max_iter * 0.015), bar.get_y() + bar.get_height()/2,
                     f'{width:.1f}',
                     va='center', ha='left', fontsize=10, fontweight='bold', color='#2c3e50')
            
        ax2.set_xlabel('Average Iterations per Timestep', fontsize=12, fontweight='bold', labelpad=10)
        ax2.set_ylabel('Configuration', fontsize=12, fontweight='bold', labelpad=10)
        ax2.set_title('Solver & Preconditioner Sweep Average Iterations', fontsize=14, fontweight='bold', pad=15)
        ax2.grid(True, linestyle='--', alpha=0.5, axis='x')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_color('#bdc3c7')
        ax2.spines['bottom'].set_color('#bdc3c7')
        plt.tight_layout()
        output_image2 = os.path.join(directory, "average_iterations.png")
        plt.savefig(output_image2, dpi=300)
        plt.close(fig2)
        print(f"Successfully generated average iterations plot: {output_image2}")
        
    except ImportError:
        print("Matplotlib is not installed. Generating text-based ASCII plot instead:")
        text_plot(data)
        print("Tip: Install matplotlib (`pip install matplotlib`) to generate high-quality PNG plots.")

if __name__ == '__main__':
    main()
