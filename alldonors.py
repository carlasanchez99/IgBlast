import os
import subprocess
import sys  # Para leer argumentos de la línea de comandos
import shutil

# Verificar que se pase la ruta como argumento
if len(sys.argv) < 2:
    print("Uso: python alldonors.py <ruta_donors>")
    sys.exit(1)

# Ruta de la carpeta principal (donors) pasada como argumento
donors_path = os.path.normpath(sys.argv[1])  # Normaliza la ruta

# Ruta de la carpeta general "Tables" dentro de "MAL67"
tables_folder = os.path.join(donors_path, "Tables")
os.makedirs(tables_folder, exist_ok=True)  # Crear la carpeta si no existe

# Lista de subcarpetas (pacientes)
pacientes = os.listdir(donors_path)

# Recorre cada subcarpeta (paciente)
for paciente in pacientes:
    paciente_path = os.path.join(donors_path, paciente)
    
    # Verifica si es una subcarpeta
    if os.path.isdir(paciente_path):
        print(f"Procesando paciente {paciente}...")

        # Construir rutas específicas para cada paciente
        input_folder = os.path.join(paciente_path, "Igblastn_output_noextend")  # Carpeta de entrada
        fasta_folder = os.path.join(paciente_path, "Sequences_fasta_aa")        # Secuencias FASTA
        igblastCregion5 = os.path.join(paciente_path, "Igblastn_output")        # Salida Igblast
        output_excel_path = os.path.join(paciente_path, f"LM_{paciente}_threshold80.xlsx")  # Salida por paciente
        general_excel_path = os.path.join(tables_folder, f"LM_{paciente}_threshold80.xlsx")  # Salida en carpeta general

        # Normalizar rutas
        input_folder = os.path.normpath(input_folder)
        fasta_folder = os.path.normpath(fasta_folder)
        igblastCregion5 = os.path.normpath(igblastCregion5)
        output_excel_path = os.path.normpath(output_excel_path)
        general_excel_path = os.path.normpath(general_excel_path)

        # Verificar que las carpetas necesarias existan
        if os.path.exists(input_folder) and os.path.exists(fasta_folder) and os.path.exists(igblastCregion5):
            print(f"Procesando archivos para el paciente {paciente}...")
            try:
                # Llama a parser_carla.py
                subprocess.run(['python', 'parser_carla.py', input_folder, fasta_folder, igblastCregion5, output_excel_path], check=True)
                
                # Copiar el archivo generado a la carpeta general
                if os.path.exists(output_excel_path):
                    os.makedirs(tables_folder, exist_ok=True)
                    shutil.copy(output_excel_path, general_excel_path)
                    print(f"Tabla de {paciente} guardada también en la carpeta general.")
                
                print(f"Paciente {paciente} procesado correctamente.")
            except subprocess.CalledProcessError as e:
                print(f"Error procesando el paciente {paciente}: {e}")
        else:
            print(f"Carpetas necesarias para el paciente {paciente} no están completas.")
