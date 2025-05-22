import sys
import numpy as np
from itertools import combinations

def find_cys_residues(pdb_file:str, bfactor:int=35, pLDDT:int=40, is_alfaphold:bool=False):
    """
    Extrae los residuos de cisteína (C, CYS) de un archivo PDB, filtrando por valores de incertidumbre 
    (B-factor, para estructuras experimentales; pLDDT para predicciones de Alphafold)

    Args:
        - pdb_file (str): Ruta al archivo PDB de la proteína
        - bfactor (int, opcional): Umbral máximo para el B-factor (valor por defecto 35)
        - pLDDT (int, opcional): Umbral mínimo para pLDDT en estructuras Alphafold (por defecto 40)
        - is_alfaphold (bool, opcional): Indica si la estructura procede de Alphafold, para filtrar con pLDDT

    Returns:
        - dict: Diccionario de residuos de CYS donde las claves son tuplas (cadena, número de residuo) y los valores
        son subdiccionarios con las coordenadas (x,y,z) de los átomos de ese residuo
        (or) 
        - None: Si no se consigue extraer ningún residuo, no se devuelve nada
    """

    cys_data = {}
    infile = open(pdb_file, 'r')

    for line in infile:
        line = line.strip('\n')

        if line.startswith('ATOM') and line[17:20].strip() == "CYS":
            atom_name = line[12:16].strip() # Nombre del átomo
            res_id = (line[21], int(line[22:26])) # Cadena, Número de residuo
            atom_x = float(line[30:38]) # Coordenada en x del átomo
            atom_y = float(line[38:46]) # Coordenada en y del átomo
            atom_z = float(line[46:54]) # Coordenada en z del átomo
            uncertainty = float(line[60:66]) # Valor del parámetro de incertidumbre

            # Criterios de filtrado de átomos para AF y Exp
            if is_alfaphold:
                if uncertainty < pLDDT:
                    continue
            else:
                if uncertainty > bfactor:
                    continue
            
            if res_id not in cys_data:
                cys_data[res_id] = {}

            cys_data[res_id][atom_name] = (atom_x, atom_y, atom_z)
    infile.close()

    return cys_data if len(cys_data) > 0 else None

def calc_distance(a, b):
    """
    Calcula la distancia entre dos puntos en tres dimensiones. La función asume
    que se le van a proporcionar objetos de tipo np.ndaarray (numpy) para describir
    los puntos a y b

    Args:
        - a (np.ndarray): Coordenadas del primer punto (x, y, z)
        - b (np.ndarray): Coordenadas del segundo punto (x, y, z)

    Returns:
        - float: Distancia entre los dos puntos
    """

    return np.linalg.norm(a - b)

def calc_dihedral(p1, p2, p3, p4):
    """
    Calcula el ángulo diedro entre los planos formados por: (CB1, SG1, SG2) y (SG1, SG2, CB2). 
    Se empleará posteriormente en la predicción de la formación de puentes disulfuro (S-S)

    Args:
        - p1, p2, p3, p4 (np.ndarray): Coordenadas de los cuatro puntos (x,y,z)

    Returns:
        - (float): Ángulo diedro en grados
    """
    
    b0 = p2 - p1 # Vector entre CB1 y SG1
    b1 = p3 - p2 # Vector entre SG1 y SG2 (enlace disulfuro)
    b2 = p4 - p3 # Vector entre SG2 y CB2 

    b1 /= np.linalg.norm(b1) # Normalizamos b1 (módulo 1, vector unitario)

    # Proyecciones ortogonales
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # Producto escalar entre los vectores proyectados
    x = np.dot(v, w)

    # Determinamos la posición relativa de los planos
    y = np.dot(np.cross(b1, v), w)

    # Devolvemos el ángulo en grados
    return np.degrees(np.arctan2(y, x))

def find_disulfide_candidates(cys_data, angle_range=(84, 96), dist_range=(1.5, 2.5)):
    """
    Busca pares de residuos de CYS que podrían formar puentes disulfuro basándose en
    la distancia entre átomos de azufre (SG1, SG2) y el ángulo diedro que forman (χSS)

    Args:
        - cys_data (dict): Diccionario con coordenadas de átomos CYS
        - angle_range (tuple, opcional): Rango permitido para el ángulo diedro en grados (por defecto, 84-96 grados angulares)
        - dist_range (tuple, opcional): Rango permitido para la distancia S-S en Ångströms (por defecto, 1.5-2.5)

    Returns:
        - list: Lista de tuplas (res1, res2, distancia, ángulo χSS) de pares candidatos
        (or)
        - None: si no se encuentran candidatos
    """

    results = []

    for (res1, atoms1), (res2, atoms2) in combinations(cys_data.items(), 2):
        if all(atom in atoms1 for atom in ('SG', 'CB')) and all(atom in atoms2 for atom in ('SG', 'CB')):
            # Extraemos las coordenadas de los átomos SG y CB de los residuos involucrados
            SG1 = np.array(atoms1['SG'])
            SG2 = np.array(atoms2['SG'])
            CB1 = np.array(atoms1['CB'])
            CB2 = np.array(atoms2['CB'])

            # Calculamos la distancia entre los átomos de azufre de las cisteínas
            dist = calc_distance(SG1, SG2)

            # Si está fuera del rango, se descarta la pareja
            if not (dist_range[0] <= dist <= dist_range[1]):
                continue
            
            # Si la distancia es aceptable, se evalúa el ángulo diedro
            chi_ss = calc_dihedral(CB1, SG1, SG2, CB2)

            # Si el ángulo también es aceptable, se guarda la combinación de residuos
            if angle_range[0] <= abs(chi_ss) <= angle_range[1]:
                results.append((res1, res2, dist, chi_ss))

    return results if len(results) > 0 else None

def main():
    """
    Programa principal.

    Organiza la ejecución del programa para detectar posibles puentes disulfuro 
    entre residuos de cisteína (CYS) en una estructura proteínica en formato PDB.
    Puede funcionar tanto con modelos experimentales como con predicciones de AlphaFold, 
    usando umbrales de calidad estructural específicos (B-factor o pLDDT).

    Uso: python3 <disulfide_find.py> <protein.pdb> [-a or --alphafold]
        - <disulfide_find.py> Nombre del programa
        - <protein.pdb> Ruta al archivo de la estructura de la proteína de interés en formato PDB
        - [--alphafold] Argumento opcional para indicar si su estructura procede de una predicción de Alphafold

    @author: Pablo López Nieto
    @date: 22-05-2025

    """
    if len(sys.argv) < 2:
        print("Uso del comando: python3 <disulfide_find.py> <protein.pdb> [--alphafold]"
              + "\n\t<disulfide_find.py> Nombre del programa"
              + "\n\t<protein.pdb> Ruta al archivo de la estructura de su proteína de interés en formato PDB"
              + "\n\t[--alphafold] Argumento opcional para indicar si su estructura procede de una predicción de Alphafold")
        sys.exit(1)

    pdb_file = sys.argv[1]
    is_alphafold = '--alphafold' in sys.argv[2:]
    
    if is_alphafold:
        print("AVISO: Modo AlphaFold activado")

    # Extraemos las cisteínas con suficiente certeza del archivo de entrada
    cys_data = find_cys_residues(pdb_file, is_alfaphold=is_alphafold)
    
    if cys_data:
        resultados = find_disulfide_candidates(cys_data)

        if resultados:
            # Reportamos los puentes disulfuro detectados (si existen)
            for res1, res2, dist, chi_ss in resultados:
                print(f"Potencial puente disulfuro entre los residuos {res1} y {res2} | Dist: {dist:.2f} Å | χSS: {chi_ss:.2f}°")
        else:
            print('No se han encontrado residuos de CYS compatibles espacialmente.')
    
    else:
        print('No se han encontrado residuos de CYS de suficiente calidad.')

if __name__ == '__main__':
    main()