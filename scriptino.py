import json
import sys

# Controllo dell'argomento
if len(sys.argv) != 2:
    print("Uso: python3 modify_json.py <file_json>")
    sys.exit(1)

file_path = sys.argv[1]

# Funzione per modificare i valori del JSON
def modify_json(file_path, modifications):
    try:
        # Leggere il file JSON esistente
        with open(file_path, "r") as f:
            data = json.load(f)

        # Apportare le modifiche
        for section, updates in modifications.items():
            if section in data:
                data[section].update(updates)
            else:
                print(f"Sezione '{section}' non trovata nel file JSON.")

        # Salvare le modifiche
        with open(file_path, "w") as f:
            json.dump(data, f, indent=4)
        print(f"Modifiche salvate correttamente in '{file_path}'.")

    except Exception as e:
        print(f"Errore durante la modifica del file JSON: {e}")
        sys.exit(1)

# Modifiche da apportare
modifications = {
    "dparams": {
        "alpha": 1.0,
        "beta": 0.8
    },
    "iparams": {
        "N": 20
    },
    "sparams": {
        "ddata": "new_densita.dat"
    },
    "bparams": {
        "save_conditions": True
    }
}

# Eseguire la modifica
modify_json(file_path, modifications)
