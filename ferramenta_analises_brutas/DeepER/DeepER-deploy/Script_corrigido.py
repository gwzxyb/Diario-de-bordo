import sys
import torch
import numpy as np
import os

# FORÇAR USO DE CPU
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

# Parâmetros
fasta = sys.argv[1]
parapath = sys.argv[2]
strand = sys.argv[3]
resultp = sys.argv[4]
batch = int(sys.argv[5])

print(f"Processando: {fasta}")
print(f"Batch size: {batch}")
print(f"Strand: {strand}")
print("Usando CPU (CUDA desabilitado)")

if __name__ == "__main__":
    # Carregar dados em janelas de 5000bp
    fa = quickloader(fasta, 5000, batch)
    
    # Fazer predição
    predict = PreDict(fa, DRBiLSTM, parapath, "hc", strand, True)
    resultados = predict.getallresults()
    
    # Salvar resultados em formato simples (cada linha é uma probabilidade)
    print("Salvando resultados...")
    with open(resultp, "w") as f:
        for posicao, probabilidade in resultados.items():
            # Escreve posição e probabilidade
            f.write(f"{posicao}\t{probabilidade}\n")
    
    print(f"Concluído! Resultados salvos em: {resultp}")
    print(f"Total de posições processadas: {len(resultados)}")
