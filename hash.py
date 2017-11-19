from datasketch import MinHash


def count_kmers(read, k):
    counts = {}
    num_kmers = len(read) - k + 1
    for i in range(num_kmers):
        kmer = read[i:i + k]
        if kmer not in counts:
            counts[kmer] = 0
        counts[kmer] += 1
    return counts
def getJaccardIndex(sequence1, sequence2, k):
    seq1_minHash, seq2_minHash = MinHash(), MinHash()
    seq1_kmers = count_kmers(sequence1, k)
    seq2_kmers = count_kmers(sequence2, k)
    seq1_keys = list(seq1_kmers.keys())
    seq2_keys = list(seq2_kmers.keys())
    for key in seq1_keys:
        seq1_minHash.update(key.encode('utf8'))
    for key in seq2_keys:
        seq2_minHash.update(key.encode('utf8'))
    return m2.jaccard(m1)




data1 = ['minhash', 'is', 'a', 'probabilistic', 'data', 'structure', 'for',
         'estimating', 'the', 'similarity', 'between', 'datasets']
data2 = ['minhash', 'is', 'a', 'probability', 'data', 'structure', 'for',
         'estimating', 'the', 'similarity', 'between', 'documents']

m1, m2 = MinHash(), MinHash()
# for d in data1:
# m1.update(d.encode('utf8'))
# for d in data2:
# m2.update(d.encode('utf8'))
# print("Estimated Jaccard for data1 and data2 is", m1.jaccard(m2))

s1 = set(data1)
s2 = set(data2)
actual_jaccard = float(len(s1.intersection(s2))) / float(len(s1.union(s2)))
print("Actual Jaccard for data1 and data2 is", actual_jaccard)

# x = count_kmers('MPMQDHKYFYLYSITNKTTEKIYVGVHKTSNLDDGYMGSGVAIKNAIKKYGIDNFYKHIIKFFESEKAMYDAEAEIVTEEFVKSKKTYNMKLGGIGGFPKHNTAGAKNGFYGKSHSRETRLKISIKSSRKRGPRGLEVKL', 2)
#y = count_kmers('MNYGQFEIESTIIATLLKQPDVLEKIRVKDYMFTNEKFKTFFNYVMDVGKIDHQEIYLKATKDKEFLDADTITKLYNSDFIGYGFFERYQQELLESYQLNKANELVTEFKQQPTNQNFNNLIDELKDLKTITNKKEDGTKKFVEEFVEELYSDSPKKQIKTGYKLMDYKIGGLEPSQLIVIAARPSVGKTGFALNMMLNIAQNGYKTSFFSLETTGTSVLKRMLSTITGIELTKIKEIRNLTPDDLTKLTNAMDKIMKLGIDISDKSNITPQDVRAQAMRHSDRQQVIFIDYLQLMDTDAKVDRRVAVEKISRDLKIIANETGAIIVLLSQLNRGVESRQDKRPMLSDMKESGGIEADASLAMLLYRDDYYNRDEDDSITGKSIVECNIAKNKDGETGIIEFEYYKKTQRFFT', 4)
# same cluster
#x = count_kmers('MSADRYLAELLSDEHDDGTMPDNVDALRQAVVGRRIVSATRGQAKINVNRWGGRDRLEGVTGLIIELDDGTKVILEDTSDCCAYTDLKSFLLAPDSVDHVIIGVGTTDGYETWHVYADMGDVLKLSVGWSCGNPFYYGYGFRIHVSRIIDGEIIPERKAIGS',4)
#y = count_kmers('MTDRYPVEETCDCAWGCENPDHPKGESRYCPLPEDDGTMPNNVAELAGHVVGRRIVSAQKEKVQTGRWYGESECFVITLDDGSRVALVDSGDCCAYTTLENFLLHPELVDHVITGVGTTGGYEKWHIYADLGDVLELDVSWSPGNPFYYGYGFDFVVLPVDE', 4)
#y = count_kmers('MTDRYPVEETCDCAWGCENPDHPKGESRYCPLPEDDGTMPNNVAELAGHVVGRRIVSAQKEKVQTGRWYGESECFVITLDDGSRVALVDSGDCCAYTTLENFLLHPELVDHVITGVGTTGGYERWHIYADLGDVLELDVSWSPGNPFYYGYGFDFVVLPVDE',4)

#same cluster
x = count_kmers('MTKRFFVEGAIAVWVPGNTPAYMIPKSFEWVSKDESSYFDQIPVVLVSKTAMDSVEYALPVLPDDLMEAIEGDAVVSDDDMEEGFDPDFAAILTDKNRSLIVVAEDEGYITRKSRLDPKTDYEVSSKAARYADELVEFDWDINVVEKESKEKTALEKFETPTEEVFIGITRRERKLKEILMNAIFDMGCNGNVEEIKYWISEMAPSSHTSEDLKDMDVDSLLQYLYDAVKQGWNDNHEVLGDGITRVYGGRHAKWMLFKKINKDRIVWFVY', 4)
y = count_kmers('MTKRFFVENAIAVWVPGNTPAYMIPKSFEWTPKDESSYFDEIPVVLVSKAAMDSVEYALPILPDDLMEAIAGDAVVSDDDMEEGFDPDFAAILTDKNRSLIVVAEDEGYITRKSRLDPETDYEVASKVARYADELTEFDWDINVVEKESREKTALEKFETPSEEVFIGITRRERKLKELLMNAIFDMGCNGNVEEIKYWISEMAPSGNTSKELKNMDVDTLLQYLYDAVKQGWNENHEVLGDGITKGYGGSHAQWNLFRKTNKDMIVRFAY', 4)
#same cluster
#not same as above
y = count_kmers('MKKNLISGSRERQKSGSRRKKSMNKIWMKIEEENKPTMYAIFTIDGTQWNVSDLVIEQEDGKDLSVKKNFQKRIKEEIGYKIVYSFEERQRMDREFQKEICKHVTYTTLYKAAIEAWKTMKPVMASLEEVDEDGE', 4)
print(x)
x_keys = list(x.keys())
y_keys = list(y.keys())
for key in x_keys:
    m1.update(key.encode('utf8'))
for key in y_keys:
    m2.update(key.encode('utf8'))

print("Estimated Jaccard for x and y is", m2.jaccard(m1))

s1 = set(x_keys)
s2 = set(y_keys)
actual_jaccard = float(len(s1.intersection(s2))) / float(len(s1.union(s2)))
print("Actual Jaccard for x and y is", actual_jaccard)
print(len(str(10**31)))
