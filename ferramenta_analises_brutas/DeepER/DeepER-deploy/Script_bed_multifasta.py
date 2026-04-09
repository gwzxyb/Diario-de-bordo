import sys
import torch
import numpy as np
import os
from Bio import SeqIO

# Forçar uso de CPU
os.environ['CUDA_VISIBLE_DEVICES'] = ''
torch.cuda.is_available = lambda: False

original_cuda = torch.Tensor.cuda
def fake_cuda(self, *args, **kwargs):
    return self
torch.Tensor.cuda = fake_cuda

sys.path.append('.')

from Model import DRBiLSTM
from lib import quickloader, PreDict

def processar_sequencia(sequencia, nome, parapath, strand, batch, cutoff, output_bed):
    """
    Processa uma única sequência e adiciona ao arquivo BED
    """
    # Criar arquivo FASTA temporário para esta sequência
    temp_fasta = f"temp_{nome}.fa"
    with open(temp_fasta, 'w') as f:
        f.write(f">{nome}\n")
        f.write(sequencia + "\n")
    
    # Carregar dados
    fa = quickloader(temp_fasta, 5000, batch)
    
    # Fazer predição
    predict = PreDict(fa, DRBiLSTM, parapath, "hc", strand, True)
    resultados = predict.getallresults()
    
    # Converter para BED
    todas_probs = []
    for janela_info in resultados.values():
        probs = janela_info[1]  # array de floats
        todas_probs.extend(probs)
    
    # Identificar regiões acima do cutoff
    regioes = []
    inicio_regiao = None
    
    for i, prob in enumerate(todas_probs):
        if prob >= cutoff:
            if inicio_regiao is None:
                inicio_regiao = i
        else:
            if inicio_regiao is not None:
                regioes.append((inicio_regiao, i))
                inicio_regiao = None
    
    if inicio_regiao is not None:
        regioes.append((inicio_regiao, len(todas_probs)))
    
    # Mesclar regiões próximas (distância ≤ 50bp)
    regioes_mescladas = []
    for inicio, fim in regioes:
        if not regioes_mescladas:
            regioes_mescladas.append([inicio, fim])
        else:
            ultimo_inicio, ultimo_fim = regioes_mescladas[-1]
            if inicio - ultimo_fim <= 50:
                regioes_mescladas[-1][1] = fim
            else:
                regioes_mescladas.append([inicio, fim])
    
    # Escrever no arquivo BED (modo append)
    with open(output_bed, 'a') as f:
        for i, (inicio, fim) in enumerate(regioes_mescladas):
            nome_rloop = f"R-loop_{nome}_{i+1}"
            f.write(f"{nome}\t{inicio}\t{fim}\t{nome_rloop}\t1000\t.\n")
    
    # Limpar arquivo temporário
    os.remove(temp_fasta)
    
    print(f"  {nome}: {len(regioes_mescladas)} R-loops encontrados")
    return len(regioes_mescladas)

def main():
    if len(sys.argv) != 6:
        print("Uso: python Script_bed_multifasta.py <multifasta> <DeepER.pkl> <strand> <output.bed> <cutoff>")
        print("Exemplo: python Script_bed_multifasta.py Hsalinarum.fa DeepER.pkl forward saida.bed 0.95")
        sys.exit(1)
    
    multifasta = sys.argv[1]
    parapath = sys.argv[2]
    strand = sys.argv[3]
    output_bed = sys.argv[4]
    cutoff = float(sys.argv[5])
    batch = 64
    
    print(f"Processando: {multifasta}")
    print(f"Cutoff: {cutoff}")
    print(f"Strand: {strand}")
    print("-" * 50)
    
    # Limpar arquivo de saída anterior
    if os.path.exists(output_bed):
        os.remove(output_bed)
    
    # Ler todas as sequências do FASTA
    total_rloops = 0
    for record in SeqIO.parse(multifasta, "fasta"):
        nome = record.id
        sequencia = str(record.seq)
        print(f"Processando {nome} ({len(sequencia)} pb)...")
        
        num = processar_sequencia(sequencia, nome, parapath, strand, batch, cutoff, output_bed)
        total_rloops += num
    
    print("-" * 50)
    print(f"Total: {total_rloops} R-loops encontrados em todas as sequências")
    print(f"Resultados salvos em: {output_bed}")

if __name__ == "__main__":
    main()
