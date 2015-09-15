import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import operator

# Input is a list of strings that are ints with possible white space
# Returns a list of ints
def strip_ints_list(data):
    data = [x.strip() for x in data]
    data = [int(x) for x in data if len(x) > 0]
    return data

# input is the connand window file from the fsgm program
# returns the number of new, core and total genes for each number of genomes
def get_genes_per_genome_data(my_file):
    in_core = False
    core = list()
    in_new_genes = False
    new_genes = list()
    in_total = False
    total = list()

    with open(my_file, 'rU') as f:
        for line in f:
            if line.startswith("core_stdv ="):
                 in_total = False
            if in_total:
                total.append(line)
            if line.startswith("total ="):
                in_new_genes = False
                in_total = True
            if in_new_genes:
                new_genes.append(line)
            if line.startswith("new ="):
                in_core = False
                in_new_genes = True              
            if in_core:
                core.append(line)
            if line.startswith("core ="):
                in_core = True
    return (core, new_genes, total)

# Input is a tuple with lists of counts for core, new and total genes, an
# output file name and the most likely value for n.  Writes a pdf plot to the
# output file
def plot_genes_per_genome((core, new_genes, total), out_file, best_n):
    num_genomes = range(1, len(new_genes) + 1)

    with PdfPages(out_file) as pdf:
        plt.plot(num_genomes, new_genes, 'ro', label='new')
        plt.plot(num_genomes, total, 'bs', label='total')
        plt.plot(num_genomes, core, 'g^', label='core')
        
        plt.axhline(best_n, color = 'k')
        plt.annotate(s='N=' + str(best_n), xy=(5, best_n - 120))
        
        plt.xlabel("Number of genomes")
        plt.ylabel('Number of genes')
        plt.title("Estimated size of distributed genome per strain sequenced")
        plt.legend(loc=0, numpoints=1)
        pdf.savefig()
        plt.close()

# Input is the likelihood file from the fsgm program
# returns a list of n values and a lis of their likelihoods
def get_lik_data(my_file):
    with open(my_file, 'rU') as f:
        data = f.read()
    split_data = [x.split('\t') for x in data.split('\n') if len(x) > 0]
    n = [int(x[0]) for x in split_data]
    lik = [float(x[1]) for x in split_data]
    return n, lik

# Graphs the n against lik and saves the results to out_file
# Returns the most likely value of n
def plot_lik_vs_n(n, lik, out_file):
    with PdfPages(out_file) as pdf:
        plt.plot(n, lik, 'ko')
        plt.xlabel("Number of genes")
        plt.ylabel('Log likelihood')
        max_index, max_lik = max(enumerate(lik), key=operator.itemgetter(1))
        best_n = n[max_index]
        label = '(' + str(best_n) + ' ,' + str(max_lik) + ')'
        plt.annotate(s=label, xy=(best_n, max_lik),
                     xytext=(best_n + 500, max_lik - 40),
                     arrowprops=dict(facecolor='red', shrink=0.05))
        plt.title("Number of genes in the pan genome")
        pdf.savefig()
        plt.close()
    return best_n
    
def main():

    fsgm_file = "/home/rachel/Documents/Finite_supragenome_model/fsgm/CommandWindow_pg70.txt"
    out_dir = "/home/rachel/Data/pg/roary_except_MB1843/results"
    nickname = 'testing'

    fsgm_lik_file = "/home/rachel/Documents/Finite_supragenome_model/fsgm/N_vs_likelihood_pg70.txt"
    n, lik = get_lik_data(fsgm_lik_file)
    out_file = out_dir + '/' + nickname + '_genes_in_pan_genome.pdf'
    best_n = plot_lik_vs_n(n, lik, out_file)

    genes_per_genome = get_genes_per_genome_data(fsgm_file)
    genes_per_genome = [strip_ints_list(x) for x in genes_per_genome]
    out_file = out_dir + '/' + nickname + '_new_genes_per_sequenced_genome.pdf'
    plot_genes_per_genome(genes_per_genome, out_file, best_n)
    


if __name__ == "__main__":
    main()
