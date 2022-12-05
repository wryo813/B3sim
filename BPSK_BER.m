clear;
clc;

N = 100;
SNR = linspace(0,20,N);
BER = zeros(1,N);

for index = 1:N
    %データレートよりもサンプリング周波数の方が十分大きくなければならない
    %データレート<搬送波周波数<サンプリング周波数である必要がある
    Fs = 200; %サンプリング周波数[Hz]
    Ts = 1/Fs;  %サンプリング周期[s]
    Fd = 10;    %データレート[Hz]
    Td = 1/Fd;  %データの周期（ナイキスト周期）[s]
    Fc = 50;   %搬送波周波数[Hz]
    A  = 1;     %振幅[V]
    N = 1000000;   %ランダム生成するビット数
    %TXD = [0 0 1 1 0 0 1 0 0 1];          %送信データ
    TXD = logical(randi([0, 1], [1, N]));   %送信データ(ランダム)
    t  = 0:Ts:(numel(TXD)-1)*Td;            %離散時間[s]
    
    %信号空間ダイアグラムへのマッピング
    BPSK = zeros(1,numel(TXD));
    BPSK(TXD==0) = A;
    BPSK(TXD==1) = -A;
    
    %インパルス列作成
    i = 1;
    impulse = zeros(1,numel(t));
    for data = BPSK
        impulse(i) = data;
        i = i + Fs/Fd;
    end
    
    %フィルタのインパルス応答を計算
    k = 5;%フィルタインパルス応答の打ち切り時間調整
    tf = -Td*k:Ts:Td*k-Ts; %フィルタインパルス応答の時間列
    h = (1/Td).*sinc(tf/Td);
    
    %インパルス列をフィルタに通過させる
    impulse_LPF = conv(impulse,h);
    impulse_LPF = impulse_LPF/max(impulse_LPF);
    delay = numel(h)/2; %フィルタの遅延量
    impulse_LPF(1:delay) = []; %先頭の遅延分を削除
    impulse_LPF = impulse_LPF(1:numel(t)); %後ろの遅延分を削除
    
    %インパルス列を搬送波に乗せる
    carrier = cos(2*pi*Fc*t); %搬送波
    s = impulse_LPF.*carrier;
    
    %noisePower = 10000; %ノイズ電力
    %r = s + normrnd(0,sqrt(noisePower),[1,numel(s)]);
    r = awgn(s,SNR(index));
    
    %同期検波
    rd = r.*carrier; %受信波と基準信号と掛ける
    
    %高調波を除去
    %フィルタのインパルス応答を計算（送信側のLPFと同じ）
    k = 50; %フィルタインパルス応答の打ち切り時間調整
    tf = -k*Td:Ts:Td*k-Ts;        %フィルタインパルス応答の時間列
    h = (1/Td).*sinc(tf/Td);
    
    %rdをフィルタに通過させる
    rd_LPF = conv(rd,h);
    delay = numel(h)/2;             %フィルタの遅延量
    rd_LPF(1:delay) = [];           %先頭の遅延分を削除
    rd_LPF = rd_LPF(1:numel(t));    %後ろの遅延分を削除
    
    %ビット判定（符号変換）
    RXD = zeros(1,numel(TXD));
    j = 1;
    for i = 1:numel(TXD)
        if  rd_LPF(j)>0
            RXD(i) = 0;
        else
            RXD(i) = 1;
        end
        j = j + Fs/Fd;
    end
    
    RXD; %受信ディジタル信号
    
    error = 0;
    for i = 1:numel(TXD)
        if RXD(i) ~= TXD(i)
            error = error + 1;
        end
    end
    
    error;
    BER(index) = error./numel(TXD);
end

plot(SNR,BER)