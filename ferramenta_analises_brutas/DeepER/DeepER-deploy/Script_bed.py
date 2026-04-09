import sys
import torch
import numpy as np
import os

# Forçar uso de CPU
os.environ['CUDA_VISIBLE_DEVICES'] = ''
torch.cuda.is_available = lambda: False

# Monkey patch para evitar .cuda()
original_cuda = torch.Tensor.cuda
def fake_cuda(self, *args, **kwargs):
    return self
torch.Tensor.cuda = fake_cuda

sys.path.append('.')

from Model import DRBiLSTM
from lib import quickloader, PreDict

def converter_para_bed(probabilidades_por_janela, cutoff, output_bed, nome_cromossomo="seq"):
    """
    Converte as probabilidades base a base em regiões BED
    """
    # Juntar todas as probabilidades em um único array contínuo
    todas_probs = []
    posicao_atual = 0
    
    for janela_info in probabilidades_por_janela:
        # Cada item é um array de probabilidades para uma janela de 5000bp
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
                # Fim da região (i é a primeira posição abaixo do cutoff)
                regioes.append((inicio_regiao, i))
                inicio_regiao = None
    
    # Última região se existir
    if inicio_regiao is not None:
        regioes.append((inicio_regiao, len(todas_probs)))
    
    # Mesclar regiões próximas (distância ≤ 50bp)
    regioes_mescladas = []
    for inicio, fim in regioes:
        if not regioes_mescladas:
            regioes_mescladas.append([inicio, fim])
        else:
            ultimo_inicio, ultimo_fim = regioes_mescladas[-1]
            if inicio - ultimo_fim <= 50:  # Mesclar se distância ≤ 50pb
                regioes_mescladas[-1][1] = fim
            else:
                regioes_mescladas.append([inicio, fim])
    
    # Escrever arquivo BED
    with open(output_bed, 'w') as f:
        for i, (inicio, fim) in enumerate(regioes_mescladas):
            # Formato BED: cromossomo inicio fim nome score strand
            nome = f"R-loop_DeepER_{i+1}"
            f.write(f"{nome_cromossomo}\t{inicio}\t{fim}\t{nome}\t1000\t.\n")
    
    print(f"Regiões encontradas: {len(regioes_mescladas)}")
    print(f"Arquivo BED salvo: {output_bed}")
    return regioes_mescladas

# Parâmetros
fasta = sys.argv[1]
parapath = sys.argv[2]
strand = sys.argv[3]
output_bed = sys.argv[4]
batch = int(sys.argv[5]) if len(sys.argv) > 5 else 64
cutoff = float(sys.argv[6]) if len(sys.argv) > 6 else 0.95

print(f"Processando: {fasta}")
print(f"Batch size: {batch}")
print(f"Cutoff: {cutoff}")
print(f"Strand: {strand}")

# Carregar dados
fa = quickloader(fasta, 5000, batch)

# Fazer predição
predict = PreDict(fa, DRBiLSTM, parapath, "hc", strand, True)
resultados = predict.getallresults()

# Converter para BED
# resultados é um dicionário: posição -> (sequencia, array_probs)
lista_resultados = list(resultados.values())
converter_para_bed(lista_resultados, cutoff, output_bed, "Hsalinarum")

print("Concluído!")
