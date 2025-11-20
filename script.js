// DOM 요소
const dnaInput = document.getElementById('dna-input');
const dnaLengthInput = document.getElementById('dna-length-input');
const decreaseLengthBtn = document.getElementById('decrease-length-btn');
const increaseLengthBtn = document.getElementById('increase-length-btn');
const randomDnaPreview = document.getElementById('random-dna-preview');
const analyzeBtn = document.getElementById('analyze-btn');
const clearBtn = document.getElementById('clear-btn');
const clearHistoryBtn = document.getElementById('clear-history-btn');
const generateRandomBtn = document.getElementById('generate-random-btn');
const errorMessage = document.getElementById('error-message');
const resultsSection = document.getElementById('results-section');
const historyList = document.getElementById('history-list');
const modeTabs = document.querySelectorAll('.mode-tab');
const manualInputSection = document.getElementById('manual-input-section');
const randomInputSection = document.getElementById('random-input-section');

// 코돈-아미노산 매핑 테이블 (3-letter code 사용)
const codonTable = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
};

const stopCodons = ['UAA', 'UAG', 'UGA'];

// DNA 서열 분석 함수
function analyzeDNA(dnaSequence) {
    // 입력 검증 및 대문자 변환
    const dna = dnaSequence.trim().toUpperCase().replace(/\s/g, '');
    
    // 유효한 DNA 염기 확인
    const validBases = /^[ATGC]+$/;
    if (!validBases.test(dna)) {
        throw new Error('DNA 서열은 A, T, G, C만 포함해야 합니다.');
    }
    
    if (dna.length === 0) {
        throw new Error('DNA 서열이 비어있습니다.');
    }
    
    // 1. DNA 길이 계산
    const dnaLength = dna.length;
    
    // 2. A/T/G/C 비율 계산
    const counts = {
        'A': (dna.match(/A/g) || []).length,
        'T': (dna.match(/T/g) || []).length,
        'G': (dna.match(/G/g) || []).length,
        'C': (dna.match(/C/g) || []).length
    };
    
    const percentages = {
        'A': (counts['A'] / dnaLength) * 100,
        'T': (counts['T'] / dnaLength) * 100,
        'G': (counts['G'] / dnaLength) * 100,
        'C': (counts['C'] / dnaLength) * 100
    };
    
    // 3. RNA 전사 결과 (T → U 변환)
    const rnaSequence = dna.replace(/T/g, 'U');
    
    // 4. 번역: 개시코돈(AUG)부터 STOP 코돈까지 protein 서열
    const translationResult = translateRNAToProtein(rnaSequence);
    
    return {
        dnaLength,
        baseCounts: counts,
        basePercentages: percentages,
        rnaSequence,
        proteinSequence: translationResult.proteinSequence,
        startCodonIndex: translationResult.startCodonIndex,
        stopCodonIndex: translationResult.stopCodonIndex,
        hasStartCodon: translationResult.hasStartCodon,
        hasStopCodon: translationResult.hasStopCodon
    };
}

// RNA를 단백질로 번역하는 함수
function translateRNAToProtein(rnaSequence) {
    let protein = [];
    let startFound = false;
    let startCodonIndex = -1;
    let stopCodonIndex = -1;
    
    // AUG(개시코돈)을 찾을 때까지 탐색
    let i = 0;
    while (i < rnaSequence.length - 2) {
        const codon = rnaSequence.substring(i, i + 3);
        
        // 개시코돈 찾기
        if (codon === 'AUG' && !startFound) {
            startFound = true;
            startCodonIndex = i;
            protein.push('Met'); // AUG는 메티오닌(Met)으로 번역
            i += 3;
            continue;
        }
        
        // 개시코돈을 찾은 후 번역 계속
        if (startFound) {
            // STOP 코돈 확인
            if (stopCodons.includes(codon)) {
                stopCodonIndex = i;
                break;
            }
            
            // 코돈을 아미노산으로 변환
            if (codonTable[codon]) {
                protein.push(codonTable[codon]);
            } else {
                // 잘못된 코돈은 'Xxx'로 표시
                protein.push('Xxx');
            }
            
            i += 3;
        } else {
            i += 1;
        }
    }
    
    // 개시코돈이 없는 경우
    if (!startFound) {
        return {
            proteinSequence: '개시코돈이 없어 아미노산 합성이 이루어지지 않습니다.',
            startCodonIndex: -1,
            stopCodonIndex: -1,
            hasStartCodon: false,
            hasStopCodon: false
        };
    }
    
    // 개시코돈은 찾았지만 단백질이 없는 경우
    if (protein.length === 0) {
        return {
            proteinSequence: '번역 가능한 서열이 없습니다.',
            startCodonIndex: startCodonIndex,
            stopCodonIndex: stopCodonIndex,
            hasStartCodon: true,
            hasStopCodon: stopCodonIndex !== -1
        };
    }
    
    // 종결코돈이 없는 경우 체크
    const hasStopCodon = stopCodonIndex !== -1;
    
    // 3-letter code를 '-'로 구분하여 반환
    let proteinSeq = protein.join('-');
    
    // 종결코돈이 없으면 경고 메시지 추가
    if (!hasStopCodon) {
        proteinSeq += '\n\n※ 종결코돈이 없어 아미노산 합성이 무한히 이어집니다.';
    }
    
    return {
        proteinSequence: proteinSeq,
        startCodonIndex: startCodonIndex,
        stopCodonIndex: stopCodonIndex,
        hasStartCodon: true,
        hasStopCodon: hasStopCodon
    };
}

// 결과 표시 함수
function displayResults(results) {
    // 1. DNA 길이
    document.getElementById('dna-length').textContent = `${results.dnaLength} bp`;
    
    // 2. A/T/G/C 비율
    const baseRatiosDiv = document.getElementById('base-ratios');
    baseRatiosDiv.innerHTML = '';
    
    ['A', 'T', 'G', 'C'].forEach(base => {
        const baseItem = document.createElement('div');
        baseItem.className = 'base-item';
        baseItem.innerHTML = `
            <div class="base-label">${base}</div>
            <div class="base-count">${results.baseCounts[base]}개</div>
            <div class="base-percentage">${results.basePercentages[base].toFixed(2)}%</div>
        `;
        baseRatiosDiv.appendChild(baseItem);
    });
    
    // 3. RNA 전사 결과 (AUG와 종결코돈 하이라이트)
    const rnaDisplay = formatRNASequence(results.rnaSequence, results.startCodonIndex, results.stopCodonIndex);
    document.getElementById('rna-sequence').innerHTML = rnaDisplay;
    
    // 4. 단백질 서열
    const proteinSeqEl = document.getElementById('protein-sequence');
    const proteinSeq = results.proteinSequence || '-';
    
    // 종결코돈 경고 메시지가 포함되어 있는지 확인
    if (proteinSeq.includes('※')) {
        const parts = proteinSeq.split('※');
        proteinSeqEl.innerHTML = `<span class="protein-sequence-main">${parts[0].trim()}</span><br><span class="protein-warning">※${parts[1]}</span>`;
    } else {
        proteinSeqEl.textContent = proteinSeq;
    }
    
    // 결과 섹션 표시
    resultsSection.style.display = 'block';
}

// 에러 메시지 표시 함수
function showError(message) {
    errorMessage.textContent = message;
    errorMessage.style.display = 'block';
    resultsSection.style.display = 'none';
}

// 에러 메시지 숨기기 함수
function hideError() {
    errorMessage.style.display = 'none';
}

// RNA 서열 포맷팅 함수 (5', 3' 표시 및 AUG, 종결코돈 하이라이트)
function formatRNASequence(rnaSequence, startCodonIndex, stopCodonIndex) {
    if (!rnaSequence) return '-';
    
    let formatted = '5\' ';
    let i = 0;
    
    while (i < rnaSequence.length) {
        // AUG(개시코돈) 위치 확인
        if (i === startCodonIndex && startCodonIndex !== -1) {
            formatted += '<span class="start-codon">AUG</span>';
            i += 3;
            continue;
        }
        
        // 종결코돈 위치 확인
        if (i === stopCodonIndex && stopCodonIndex !== -1) {
            const stopCodon = rnaSequence.substring(i, i + 3);
            formatted += `<span class="stop-codon">${stopCodon}</span>`;
            i += 3;
            continue;
        }
        
        // 일반 염기 출력
        formatted += rnaSequence[i];
        i++;
    }
    
    formatted += ' 3\'';
    
    return formatted;
}

// localStorage에서 히스토리 가져오기
function getHistory() {
    const history = localStorage.getItem('pcrHelperHistory');
    return history ? JSON.parse(history) : [];
}

// localStorage에 히스토리 저장하기
function saveHistory(dnaSequence) {
    let history = getHistory();
    
    // 중복 제거 (이미 존재하면 제거)
    history = history.filter(item => item !== dnaSequence);
    
    // 맨 앞에 추가 (최신 것이 위에)
    history.unshift(dnaSequence);
    
    // 최대 10개만 저장
    if (history.length > 10) {
        history = history.slice(0, 10);
    }
    
    localStorage.setItem('pcrHelperHistory', JSON.stringify(history));
}

// 히스토리 표시 함수
function displayHistory() {
    const history = getHistory();
    historyList.innerHTML = '';
    
    if (history.length === 0) {
        historyList.innerHTML = '<p style="color: #999; padding: 20px; text-align: center;">저장된 기록이 없습니다.</p>';
        return;
    }
    
    history.forEach((item, index) => {
        const historyItem = document.createElement('div');
        historyItem.className = 'history-item';
        historyItem.innerHTML = `
            <div class="history-item-text">${item.substring(0, 50)}${item.length > 50 ? '...' : ''}</div>
        `;
        
        // 클릭 시 해당 서열을 입력란에 넣고 분석
        historyItem.addEventListener('click', () => {
            // 수동 입력 모드로 전환
            switchMode('manual');
            dnaInput.value = item;
            analyzeDNASequence();
        });
        
        historyList.appendChild(historyItem);
    });
}

// 랜덤 DNA 생성 함수
function generateRandomDNA(length) {
    const bases = ['A', 'T', 'G', 'C'];
    let randomDNA = '';
    
    for (let i = 0; i < length; i++) {
        randomDNA += bases[Math.floor(Math.random() * bases.length)];
    }
    
    return randomDNA;
}

// 탭 전환 기능
function switchMode(mode) {
    modeTabs.forEach(tab => {
        if (tab.dataset.mode === mode) {
            tab.classList.add('active');
        } else {
            tab.classList.remove('active');
        }
    });
    
    if (mode === 'manual') {
        manualInputSection.classList.add('active');
        randomInputSection.classList.remove('active');
        
        // DNA 입력칸으로 스크롤
        setTimeout(() => {
            dnaInput.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
            dnaInput.focus();
        }, 100);
    } else {
        manualInputSection.classList.remove('active');
        randomInputSection.classList.add('active');
    }
}

// 랜덤 DNA 생성 및 미리보기에 표시
function generateAndSetRandomDNA() {
    const length = parseInt(dnaLengthInput.value);
    
    if (!length || length < 1 || length > 10000) {
        showError('DNA 길이는 1부터 10000 사이의 숫자여야 합니다.');
        return;
    }
    
    const randomDNA = generateRandomDNA(length);
    randomDnaPreview.value = randomDNA;
    
    // 에러 메시지 숨기기
    hideError();
}

// 분석 실행 함수
function analyzeDNASequence() {
    hideError();
    
    let dnaSequence = '';
    const currentMode = document.querySelector('.mode-tab.active').dataset.mode;
    
    // 현재 모드에 따라 DNA 서열 가져오기
    if (currentMode === 'manual') {
        dnaSequence = dnaInput.value.trim();
        
        if (!dnaSequence) {
            showError('DNA 서열을 입력해주세요.');
            return;
        }
    } else {
        // 랜덤 모드
        dnaSequence = randomDnaPreview.value.trim();
        
        if (!dnaSequence) {
            showError('먼저 랜덤 DNA를 생성해주세요.');
            return;
        }
    }
    
    try {
        const results = analyzeDNA(dnaSequence);
        displayResults(results);
        
        // 히스토리에 저장
        saveHistory(dnaSequence);
        displayHistory();
    } catch (error) {
        showError(error.message);
    }
}

// 초기화 함수
function clearAll() {
    dnaInput.value = '';
    randomDnaPreview.value = '';
    dnaLengthInput.value = '100';
    hideError();
    resultsSection.style.display = 'none';
    
    // 현재 활성 탭에 따라 포커스 설정
    const currentMode = document.querySelector('.mode-tab.active').dataset.mode;
    if (currentMode === 'manual') {
        dnaInput.focus();
    } else {
        dnaLengthInput.focus();
    }
}

// 히스토리 삭제 함수
function clearHistory() {
    if (confirm('모든 입력 기록을 삭제하시겠습니까?')) {
        localStorage.removeItem('pcrHelperHistory');
        displayHistory();
    }
}

// DNA 길이 증가/감소 함수
function adjustDNALength(amount) {
    const currentValue = parseInt(dnaLengthInput.value) || 100;
    const newValue = Math.max(1, Math.min(10000, currentValue + amount));
    dnaLengthInput.value = newValue;
    
    // 입력 이벤트 발생 (다른 곳에서 값을 사용할 수 있도록)
    dnaLengthInput.dispatchEvent(new Event('input'));
}

// 이벤트 리스너
analyzeBtn.addEventListener('click', analyzeDNASequence);
clearBtn.addEventListener('click', clearAll);
clearHistoryBtn.addEventListener('click', clearHistory);
generateRandomBtn.addEventListener('click', generateAndSetRandomDNA);

// DNA 길이 조절 버튼
decreaseLengthBtn.addEventListener('click', () => adjustDNALength(-50));
increaseLengthBtn.addEventListener('click', () => adjustDNALength(50));

// 탭 전환 이벤트 리스너
modeTabs.forEach(tab => {
    tab.addEventListener('click', () => {
        switchMode(tab.dataset.mode);
    });
});

// Enter 키로 분석 (Shift+Enter는 줄바꿈)
dnaInput.addEventListener('keydown', (e) => {
    if (e.key === 'Enter' && !e.shiftKey) {
        e.preventDefault();
        analyzeDNASequence();
    }
});

// 페이지 로드 시 히스토리 표시
displayHistory();

