def analyze_dna(dna_sequence):
    """
    DNA 서열을 분석하여 다음 정보를 반환:
    1. DNA 길이
    2. A/T/G/C 비율
    3. RNA 전사 결과
    4. 번역된 단백질 서열
    """
    # 입력 검증 및 대문자 변환
    dna = dna_sequence.strip().upper().replace(' ', '')
    
    # 유효한 DNA 염기 확인
    valid_bases = {'A', 'T', 'G', 'C'}
    if not all(base in valid_bases for base in dna):
        raise ValueError("DNA 서열은 A, T, G, C만 포함해야 합니다.")
    
    if len(dna) == 0:
        raise ValueError("DNA 서열이 비어있습니다.")
    
    # 1. DNA 길이 계산
    dna_length = len(dna)
    
    # 2. A/T/G/C 비율 계산
    counts = {
        'A': dna.count('A'),
        'T': dna.count('T'),
        'G': dna.count('G'),
        'C': dna.count('C')
    }
    
    percentages = {
        'A': (counts['A'] / dna_length) * 100,
        'T': (counts['T'] / dna_length) * 100,
        'G': (counts['G'] / dna_length) * 100,
        'C': (counts['C'] / dna_length) * 100
    }
    
    # 3. RNA 전사 결과 (T → U 변환)
    rna_sequence = dna.replace('T', 'U')
    
    # 4. 번역: 개시코돈(AUG)부터 STOP(UAA, UAG, UGA)까지 protein 서열
    protein_sequence = translate_rna_to_protein(rna_sequence)
    
    return {
        'dna_length': dna_length,
        'base_counts': counts,
        'base_percentages': percentages,
        'rna_sequence': rna_sequence,
        'protein_sequence': protein_sequence
    }


def translate_rna_to_protein(rna_sequence):
    """
    RNA 서열을 단백질 서열로 번역
    개시코돈(AUG)부터 STOP 코돈(UAA, UAG, UGA)까지 변환
    """
    # 코돈-아미노산 매핑 테이블
    codon_table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    stop_codons = {'UAA', 'UAG', 'UGA'}
    protein = []
    start_found = False
    
    # AUG(개시코돈)을 찾을 때까지 탐색
    i = 0
    while i < len(rna_sequence) - 2:
        codon = rna_sequence[i:i+3]
        
        # 개시코돈 찾기
        if codon == 'AUG' and not start_found:
            start_found = True
            protein.append('M')  # AUG는 메티오닌(M)으로 번역
            i += 3
            continue
        
        # 개시코돈을 찾은 후 번역 계속
        if start_found:
            # STOP 코돈 확인
            if codon in stop_codons:
                break
            
            # 코돈을 아미노산으로 변환
            if codon in codon_table:
                amino_acid = codon_table[codon]
                protein.append(amino_acid)
            else:
                # 잘못된 코돈은 'X'로 표시
                protein.append('X')
            
            i += 3
        else:
            i += 1
    
    if not start_found:
        return "개시코돈(AUG)을 찾을 수 없습니다."
    
    if len(protein) == 0:
        return "번역 가능한 서열이 없습니다."
    
    return ''.join(protein)


def print_analysis_results(results):
    """
    분석 결과를 보기 좋게 출력
    """
    print("=" * 60)
    print("DNA 서열 분석 결과")
    print("=" * 60)
    
    # 1. DNA 길이
    print(f"\n1. DNA 길이: {results['dna_length']} bp")
    
    # 2. A/T/G/C 비율
    print("\n2. A/T/G/C 비율:")
    counts = results['base_counts']
    percentages = results['base_percentages']
    
    for base in ['A', 'T', 'G', 'C']:
        print(f"   {base}: {counts[base]}개 ({percentages[base]:.2f}%)")
    
    # 3. RNA 전사 결과
    print(f"\n3. RNA 전사 결과:")
    print(f"   {results['rna_sequence']}")
    
    # 4. 번역된 단백질 서열
    print(f"\n4. 번역된 단백질 서열:")
    print(f"   {results['protein_sequence']}")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    # 예제 실행
    print("DNA 서열을 입력하세요:")
    print("(예: ATGCGATCGATCGATCGATCGATCGATCG)")
    
    try:
        dna_input = input("\nDNA 서열: ")
        results = analyze_dna(dna_input)
        print_analysis_results(results)
    except ValueError as e:
        print(f"\n오류: {e}")
    except Exception as e:
        print(f"\n예상치 못한 오류가 발생했습니다: {e}")

