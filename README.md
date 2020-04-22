# Bioinformatics-tool
Bioinformatická aplikace
Tato bioiformatická aplikace obsahuje základní metody z oblasti sekvenační strukturální bioinformatiky, jako je analýza souborů nebo porovnávání sekvencí a struktur. Obsahuje následující: 

* Parsování FASTA souboru
* Výpočet podobností dvou sekvencí na základě Hammingovy vzdálenosti
* Sekvenční alignment vytvoření na základě Editační vzdálenosti
* Parsování PDB souboru
* Zpracování Multiple sequence alignmentu
* Určení konzervovaných oblastí v Multiple sequence alignmentu
* Určování superpozice

## Download
Pomocí git clone http://github.com/korenal/    si stáhnete tento projekt

## Požadavky
Tato aplikace jde spustit pomocí jakéhokoliv programovacího prostředí podporující Python 3.7

## Použití aplikace
Po spuštění aplikace se Vám zobrazí následující GUI uživatelské rozhraní:
![Uživatelské rozhraní](https://github.com/korenal/Bioinformatics-tool/blob/master/Bioinformatics_tool_user.png)

kde máte výběr možností, které tato aplikace umí, vyberte si jednu a pokračujte dál pomocí stisknutí tlačítka "OK".

### 1. FASTA parser
* Zobrazí se dialogové okno a požádá Vás o nahrání souboru, který je ve FASTA formatu (př. fasta_file.txt)
* Program vypíše všechny sekvence, které se v nahraném souboru nachází
* vyberte si jednu z nabízených sekvencí
* Program vypíše informace o dané sekvenci (popis, sekvenci, délku sekvence) a můžete si nechat vytisknout i podsekvenci

### 2. Hamming distance
* zde program nabídne možnost buď zadání sekvencí pomocí FASTA souboru (např. fasta_file.txt), nebo zadání sekvencí ručně
* je nutné, aby sekvence měly stejnou délku, jinak program vypíše chybovou hlášku (The two sequences don´t have a same length), pokud mají dvě sekvence stejnou délku, tak vypíše jejich Hammingovu vzdálenost

### 3. Sequence alignment using edit distance
* Zde napište dvě sekvence a program spočítá jejich editační vzdálenost a vypíše všechny možné zarovnání těchto sekvencí

### 4. PDB parser
* Zobrazí se dialogové okno a požádá Vás o nahrání souboru, který je v PDB formatu (př. pdb_file.txt)
* zde máte na výběr:
1. Výpis modelu
1. Výpis struktury modelu
1. Výpis reziduí ve struktuře modelu
1. Výpis atomů v reziduích
1. Výpis informací o 
