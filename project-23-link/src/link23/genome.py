"""Approximate genetic map lengths (centimorgans) for 22 autosomes + X.

These are rounded, sex-averaged orders of magnitude, good enough for a
teaching simulation. Replace them with a real map (e.g. HapMap / deCODE)
before drawing quantitative conclusions.
"""

CHROM_NAMES = [str(i) for i in range(1, 23)] + ["X"]

CHROM_LENGTHS_CM = {
    "1": 280, "2": 260, "3": 220, "4": 210, "5": 200, "6": 190,
    "7": 180, "8": 165, "9": 155, "10": 170, "11": 150, "12": 160,
    "13": 120, "14": 110, "15": 110, "16": 120, "17": 120, "18": 115,
    "19": 95, "20": 95, "21": 60, "22": 70, "X": 180,
}
