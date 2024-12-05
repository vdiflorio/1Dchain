import tarfile
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from scipy.optimize import curve_fit
import matplotlib.ticker as mticker

import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Process tar.gz archive with specified N value.")
parser.add_argument("N", type=int, help="Value of N to process")
args = parser.parse_args()

# Get N value from the command line argument
N = args.N

# Construct the archive path based on N
archive_paths = [
    f"./script/compressed_archive_N_{N}.tar.gz"
]



# Funzione di fit (es. lineare)
def linear_fit(x, a, b):
    return a * x + b

# Liste per i dati del fit
y_final_values = []
divisor_values = []

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
                
                x_normalized = x / divisor*float(N)
                y_normalized = y / divisor*float(N)
                print(np.mean(x_normalized[-1:-20000:-1]))
                y_final_values.append(y[-1])
                divisor_values.append(divisor/float(N))


                # Create DataFrames
                df = pd.DataFrame(x_normalized, columns=["x_normalized"])
                
                integral_df = pd.DataFrame(y_normalized, columns=["y_normalized"])

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
                   # Fit dei dati raccolti
    divisor_values = np.array(divisor_values)
    y_final_values = np.array(y_final_values)

    try:
    # Fit lineare
        popt, pcov = curve_fit(linear_fit, divisor_values, y_final_values)
        a, b = popt
        perr = np.sqrt(np.diag(pcov))  # Incertezze sui parametri

        print(f"Fit parameters: a={a} ± {perr[0]}, b={b} ± {perr[1]}")

        # Genera i dati del fit per il plot
        x_fit = np.linspace(min(divisor_values), max(divisor_values), 500)
        y_fit = linear_fit(x_fit, a, b)

        # Plot dei risultati del fit
        plt.figure(figsize=(8, 5))
        plt.errorbar(divisor_values, y_final_values, fmt='o', label="Data", color="blue")
        plt.plot(x_fit, y_fit, label=f"Fit: y = ({a:.2e} ± {perr[0]:.2e})x + ({b:.2e} ± {perr[1]:.2e})", color="red")

        # Configurazione asse y in notazione scientifica
        plt.gca().yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

        # Etichette e titolo
        plt.xlabel("DT/DX")
        plt.ylabel("Flux")
        plt.title("Fit of y[-1] vs Divisor with Uncertainties")
        plt.legend()
        plt.grid()

        # Salva il plot
        fit_plot_path = os.path.join(output_dir, f"N_{N}_fit_df.pdf")
        plt.savefig(fit_plot_path)
        plt.close()

        print(f"Saved fit plot to {fit_plot_path}")
    except Exception as e:
        print(f"Error during fitting: {e}")


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

 