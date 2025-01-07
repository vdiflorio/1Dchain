import json
import sys

# Controllo degli argomenti
if len(sys.argv) != 4:
    print("Uso: python3 modify_n_tr.py <file_json> <nuovo_N> <nuovo_grad>")
    sys.exit(1)

file_path = sys.argv[1]
new_n = int(sys.argv[2])
new_grad = float(sys.argv[3])
new_tr = 1 + new_grad*new_n

# Funzione per modificare `N` e `Tr`
def modify_n_and_tr(file_path, new_n, new_tr):
    try:
        # Leggere il file JSON esistente
        with open(file_path, "r") as f:
            data = json.load(f)

        # Modificare i valori di `N` e `Tr`
        if "iparams" in data and "N" in data["iparams"]:
            data["iparams"]["N"] = new_n
        else:
            print("Parametro 'N' non trovato in 'iparams'.")

        if "dparams" in data and "Tr" in data["dparams"]:
            data["dparams"]["Tr"] = new_tr
        else:
            print("Parametro 'Tr' non trovato in 'dparams'.")

        # Salvare il file JSON aggiornato
        with open(file_path, "w") as f:
            json.dump(data, f, indent=4)
        print(f"File JSON modificato correttamente: N = {new_n}, Tr = {new_tr}")

    except Exception as e:
        print(f"Errore durante la modifica del file JSON: {e}")
        sys.exit(1)

# Eseguire la modifica
modify_n_and_tr(file_path, new_n, new_tr)
