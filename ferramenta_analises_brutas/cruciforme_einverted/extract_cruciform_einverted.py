#!/usr/bin/env python3
import re
import sys

def extract_cruciforms(arquivo_in, arquivo_out, min_arm=10, max_arm=100, max_gap=3):
    """
    Extrai cruciformes do output do einverted e salva em BED
    """
    
    with open(arquivo_in, 'r') as f:
        conteudo = f.read()
    
    # Dividir por blocos (cada repeat é separado por linhas em branco)
    blocos = conteudo.split('\n\n')
    
    resultados = []
    
    for bloco in blocos:
        if not bloco.strip():
            continue
        
        linhas = bloco.strip().split('\n')
        
        # Primeira linha: cabeçalho
        cabecalho = linhas[0]
        
        # Extrair cromossomo e score
        match_cab = re.search(r'^(\S+): Score (\d+)', cabecalho)
        if not match_cab:
            continue
        
        cromossomo = match_cab.group(1)
        score = int(match_cab.group(2))
        
        # Extrair gaps
        match_gap = re.search(r'(\d+) gaps', cabecalho)
        gaps = int(match_gap.group(1)) if match_gap else 0
        
        # Procurar linhas com sequências (padrão: números, letras, números)
        sequencias = []
        for linha in linhas:
            # Padrão: espaços, número, espaço, letras (com ou sem hífen), espaço, número
            match_seq = re.search(r'\s+(\d+)\s+([A-Za-z-]+)\s+(\d+)', linha)
            if match_seq:
                pos1 = int(match_seq.group(1))
                seq = match_seq.group(2).replace('-', '')  # remove hífens
                pos2 = int(match_seq.group(3))
                sequencias.append((pos1, seq, pos2))
        
        # Precisa de pelo menos 2 sequências (forward e reverse)
        if len(sequencias) >= 2:
            # Pega as duas primeiras
            pos1, seq1, pos2 = sequencias[0]
            pos3, seq2, pos4 = sequencias[1]
            
            # Comprimento do braço (sem gaps)
            arm_len = len(seq1)
            
            # Calcular tamanho do espaçador (gap)
            # A posição final da primeira sequência e início da segunda
            if pos3 > pos2:
                gap_len = pos3 - pos2 - 1
            else:
                gap_len = pos1 - pos4 - 1
            
            # Aplicar filtros
            if min_arm <= arm_len <= max_arm and gap_len <= max_gap:
                start = min(pos1, pos2, pos3, pos4)
                end = max(pos1, pos2, pos3, pos4)
                
                resultados.append({
                    'cromossomo': cromossomo,
                    'start': start,
                    'end': end,
                    'score': score,
                    'arm_len': arm_len,
                    'gap_len': gap_len,
                    'gaps': gaps
                })
    
    # Escrever arquivo BED
    with open(arquivo_out, 'w') as f:
        for i, r in enumerate(resultados, 1):
            nome = f"Cruciform_{i}"
            f.write(f"{r['cromossomo']}\t{r['start']}\t{r['end']}\t{nome}\t{r['score']}\t.\n")
    
    # Estatísticas
    print(f"\n{'='*60}")
    print(f"RESUMO DA EXTRAÇÃO DE CRUCIFORMES")
    print(f"{'='*60}")
    print(f"Total de blocos analisados: {len(blocos)}")
    print(f"Total de cruciformes encontrados: {len(resultados)}")
    print(f"Arquivo BED salvo: {arquivo_out}")
    
    if resultados:
        arm_lens = [r['arm_len'] for r in resultados]
        gap_lens = [r['gap_len'] for r in resultados]
        scores = [r['score'] for r in resultados]
        
        print(f"\n{'='*60}")
        print(f"ESTATÍSTICAS DOS CRUCIFORMES")
        print(f"{'='*60}")
        print(f"Comprimento do braço (pb):")
        print(f"  Mínimo: {min(arm_lens)}")
        print(f"  Máximo: {max(arm_lens)}")
        print(f"  Médio: {sum(arm_lens)/len(arm_lens):.1f}")
        
        print(f"\nTamanho do espaçador/gap (pb):")
        print(f"  Mínimo: {min(gap_lens)}")
        print(f"  Máximo: {max(gap_lens)}")
        print(f"  Médio: {sum(gap_lens)/len(gap_lens):.1f}")
        
        print(f"\nScore:")
        print(f"  Mínimo: {min(scores)}")
        print(f"  Máximo: {max(scores)}")
        print(f"  Médio: {sum(scores)/len(scores):.1f}")
        
        # Mostrar primeiros resultados
        print(f"\n{'='*60}")
        print(f"PRIMEIROS 10 CRUCIFORMES")
        print(f"{'='*60}")
        for i, r in enumerate(resultados[:10], 1):
            print(f"{i:3d}. {r['cromossomo']}:{r['start']}-{r['end']} | braço={r['arm_len']}pb | gap={r['gap_len']}pb | score={r['score']}")
    
    return resultados

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Uso: python extract_cruciform_einverted.py <arquivo_einverted.out> <saida.bed>")
        print("Exemplo: python extract_cruciform_einverted.py Hsalinarum_einverted.out cruciformes.bed")
        sys.exit(1)
    
    arquivo_in = sys.argv[1]
    arquivo_out = sys.argv[2]
    
    # Parâmetros opcionais
    min_arm = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    max_arm = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    max_gap = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    
    extract_cruciforms(arquivo_in, arquivo_out, min_arm, max_arm, max_gap)
