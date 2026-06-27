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
        
    # Sort by CPU count
    data.sort(key=lambda x: int(re.match(r'^(\d+)', x[1]).group(1)) if re.match(r'^(\d+)', x[1]) else 0)
        
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
        print("Usage: plot_scaling.py <results_directory>")
        sys.exit(1)
        
    directory = sys.argv[1]
    results_file = os.path.join(directory, "parallel_scalability.txt")
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
        
    # Extract CPU counts and verify format
    cpu_counts = []
    for job_id, config, time_s, avg_iter in data:
        match = re.match(r'^(\d+)', config)
        if match:
            cpu_counts.append(int(match.group(1)))
        else:
            print(f"Warning: Configuration '{config}' does not match CPU count pattern. Defaulting CPU count to 1.")
            cpu_counts.append(1)
            
    # Attempt to plot with matplotlib
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        # Modern styling
        plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
        
        # Sort CPU configurations numerically
        sorted_data = sorted(zip(cpu_counts, data), key=lambda x: x[0])
        sorted_cpus = [x[0] for x in sorted_data]
        sorted_times = [x[1][2] for x in sorted_data]
        sorted_iters = [x[1][3] for x in sorted_data]
        
        # Setup ticks and labels (CPU 0 runs sequentially)
        x_ticks = sorted_cpus
        x_labels = [str(x) if x > 0 else '0 (No MPI)' for x in x_ticks]
        
        # Plot 1: Execution Time (Line Plot)
        fig1, ax1 = plt.subplots(figsize=(8, 5))
        ax1.plot(sorted_cpus, sorted_times, marker='o', markersize=6, linewidth=2, color='#2b5c8f', label='Execution Time')
        
        # Add labels above points
        max_time = max(sorted_times) if sorted_times else 1.0
        for x, y in zip(sorted_cpus, sorted_times):
            ax1.text(x, y + (max_time * 0.03), f'{y:.2f}m', 
                     ha='center', va='bottom', fontsize=9, fontweight='bold', color='#2c3e50')
        
        ax1.set_xlabel('Number of MPI proccesses', fontsize=12, fontweight='bold', labelpad=10)
        ax1.set_ylabel('Execution Time (minutes)', fontsize=12, fontweight='bold', labelpad=10)
        ax1.set_title('Scaling on parallel infrastructure', fontsize=14, fontweight='bold', pad=15)
        ax1.grid(True, linestyle='--', alpha=0.5)
        ax1.set_xticks(x_ticks)
        ax1.set_xticklabels(x_labels)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_color('#bdc3c7')
        ax1.spines['bottom'].set_color('#bdc3c7')
        plt.tight_layout()
        
        output_image1 = os.path.join(directory, "execution_times.png")
        plt.savefig(output_image1, dpi=300)
        plt.close(fig1)
        print(f"Successfully generated execution times line plot: {output_image1}")
        
        # Plot 2: Average Iterations (Line Plot)
        fig2, ax2 = plt.subplots(figsize=(8, 5))
        ax2.plot(sorted_cpus, sorted_iters, marker='s', markersize=6, linewidth=2, color='#e74c3c', label='Average Iterations')
        
        # Add labels above points
        max_iter = max(sorted_iters) if sorted_iters else 1.0
        for x, y in zip(sorted_cpus, sorted_iters):
            ax2.text(x, y + (max_iter * 0.03 if max_iter > 0 else 0.2), f'{y:.1f}', 
                     ha='center', va='bottom', fontsize=9, fontweight='bold', color='#2c3e50')
        
        ax2.set_xlabel('Number of MPI proccesses', fontsize=12, fontweight='bold', labelpad=10)
        ax2.set_ylabel('Average Iterations per Timestep', fontsize=12, fontweight='bold', labelpad=10)
        ax2.set_title('Average Iterations based on MPI proccesses', fontsize=14, fontweight='bold', pad=15)
        ax2.grid(True, linestyle='--', alpha=0.5)
        ax2.set_xticks(x_ticks)
        ax2.set_xticklabels(x_labels)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_color('#bdc3c7')
        ax2.spines['bottom'].set_color('#bdc3c7')
        plt.tight_layout()
        
        output_image2 = os.path.join(directory, "average_iterations.png")
        plt.savefig(output_image2, dpi=300)
        plt.close(fig2)
        print(f"Successfully generated average iterations line plot: {output_image2}")
            
    except ImportError:
        print("Matplotlib is not installed. Generating text-based ASCII plot instead:")
        text_plot(data)
        print("Tip: Install matplotlib (`pip install matplotlib`) to generate high-quality PNG plots.")

if __name__ == '__main__':
    main()
