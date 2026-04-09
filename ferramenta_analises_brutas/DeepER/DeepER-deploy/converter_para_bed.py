import sys

def converter_prob_para_bed(arquivo_prob, cutoff, arquivo_bed):
    """
    Converte arquivo de probabilidades (posição \t probabilidade) para BED
    Agrupa posições consecutivas acima do cutoff em regiões
    """
    with open(arquivo_prob, 'r') as f:
        linhas = f.readlines()
    
    regioes = []
    inicio_atual = None
    fim_atual = None
    
    for i, linha in enumerate(linhas):
        if not linha.strip():
            continue
        
        # Assumindo formato: posição \t probabilidade
        partes = linha.strip().split('\t')
        if len(partes) >= 2:
            try:
                posicao = int(partes[0])
                prob = float(partes[1])
            except:
                continue
            
            if prob >= cutoff:
                if inicio_atual is None:
                    inicio_atual = posicao
                    fim_atual = posicao + 1
                else:
                    fim_atual = posicao + 1
            else:
                if inicio_atual is not None:
                    # Fim da região
                    regioes.append((inicio_atual, fim_atual))
                    inicio_atual = None
                    fim_atual = None
    
    # Adicionar última região se existir
    if inicio_atual is not None:
        regioes.append((inicio_atual, fim_atual))
    
    # Escrever arquivo BED
    with open(arquivo_bed, 'w') as f:
        for i, (inicio, fim) in enumerate(regioes):
            # Formato BED: cromossomo inicio fim nome score strand
            f.write(f"seq1\t{inicio}\t{fim}\tR-loop_{i+1}\t1000\t.\n")
    
    print(f"Encontradas {len(regioes)} regiões com probabilidade >= {cutoff}")
    print(f"Arquivo BED salvo em: {arquivo_bed}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python converter_para_bed.py <arquivo_prob.txt> <cutoff> <saida.bed>")
        sys.exit(1)
    
    converter_prob_para_bed(sys.argv[1], float(sys.argv[2]), sys.argv[3])
