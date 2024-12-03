import tarfile
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# List of archive paths
archive_paths = [
    "./single_data/compressed_archive_N_30.tar.gz"
]

# Directory principale per salvare i PDF
main_output_dir = "./plots"

for archive_path in archive_paths:
    # Estrai il valore di N dal percorso dell'archivio
    match = re.search(r'compressed_archive_N_(\d+)', archive_path)
    if match:
        N = match.group(1)
        output_dir = os.path.join(main_output_dir, f"N_{N}")  # Crea una cartella basata su N
    else:
        N = "Unknown"
        output_dir = os.path.join(main_output_dir, "Unknown_N")  # Cartella generica se N non trovato

    os.makedirs(output_dir, exist_ok=True)  # Crea la cartella per questo valore di N

    # Liste per i dati dei plot finali
    all_df_data = []
    all_integral_df_data = []

    with tarfile.open(archive_path, 'r:gz') as archive:
        # List all .dat files in the archive
        files = [member.name for member in archive.getmembers() if member.name.endswith('.dat')]

        for file_name in files:
            with archive.extractfile(file_name) as file:
                # Skip "Tr_1.dat"
                if "Tr_1.dat" in file_name:
                    continue

                # Extract Tr value
                match = re.search(r'Tr_([0-9]+(?:\.[0-9]+)?)', file_name)
                if match:
                    try:
                        Tr = float(match.group(1))  # Converte Tr in float
                    except ValueError:
                        print(f"Error: Couldn't convert Tr value '{match.group(1)}' to float in filename {file_name}, skipping...")
                        continue
                else:
                    print(f"Warning: Couldn't extract Tr from filename {file_name}, skipping...")
                    continue

                # Avoid division by zero
                if Tr <= 1:
                    print(f"Warning: Skipping file {file_name} due to Tr={Tr} <= 1")
                    continue

                # Load the data
                data = np.loadtxt(file)
                if data.shape[1] != 2:
                    print(f"Warning: {file_name} doesn't have exactly two columns, skipping...")
                    continue

                # Normalize data
                x = data[:, 0]
                y = data[:, 1]
                divisor = Tr - 1
                x_normalized = x / divisor
                y_normalized = y / divisor

                # Create DataFrames
                df = pd.DataFrame(x_normalized[:50000], columns=["x_normalized"])
                integral_df = pd.DataFrame(y_normalized[:50000], columns=["y_normalized"])

                # Salvare i dati per i plot finali
                all_df_data.append((file_name, df))
                all_integral_df_data.append((file_name, integral_df))

                # Plot df
                plt.figure(figsize=(8, 5))
                plt.plot(df, label="Normalized x")
                plt.xlabel("Index")
                plt.ylabel("Normalized x")
                plt.title(f"{file_name} - Normalized x")
                plt.legend()
                plt.grid()
                pdf_path_df = os.path.join(output_dir, f"N_{N}_Tr_{Tr}_df.pdf")  # Nominare il file come richiesto
                plt.savefig(pdf_path_df)
                plt.close()

                # Plot integral_df
                plt.figure(figsize=(8, 5))
                plt.plot(integral_df, label="Normalized y")
                plt.xlabel("Index")
                plt.ylabel("Normalized y")
                plt.title(f"{file_name} - Normalized y")
                plt.legend()
                plt.grid()
                pdf_path_integral = os.path.join(output_dir, f"N_{N}_Tr_{Tr}_integral_df.pdf")  # Nominare il file come richiesto
                plt.savefig(pdf_path_integral)
                plt.close()

                print(f"Saved plots for {file_name} to {pdf_path_df} and {pdf_path_integral}")

    # Creazione dei plot finali
    # Plot con tutte le df
    plt.figure(figsize=(10, 6))
    for file_name, df in all_df_data:
        plt.plot(df, label=file_name)
    plt.xlabel("Index")
    plt.ylabel("Normalized x")
    plt.title(f"All Normalized x (df) - N={N}")
    plt.legend(fontsize=8)
    plt.grid()
    final_df_path = os.path.join(output_dir, f"all_N_{N}_df.pdf")
    plt.savefig(final_df_path)
    plt.close()
    print(f"Saved combined plot of all df to {final_df_path}")

    # Plot con tutte le integral_df
    plt.figure(figsize=(10, 6))
    for file_name, integral_df in all_integral_df_data:
        plt.plot(integral_df, label=file_name)
    plt.xlabel("Index")
    plt.ylabel("Normalized y")
    plt.title(f"All Normalized y (integral_df) - N={N}")
    plt.legend(fontsize=8)
    plt.grid()
    final_integral_df_path = os.path.join(output_dir, f"all_N_{N}_integral_df.pdf")
    plt.savefig(final_integral_df_path)
    plt.close()
    print(f"Saved combined plot of all integral_df to {final_integral_df_path}")
