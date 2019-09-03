%% Os c�digos abaixo foram gerados tendo como refer�ncia
%% o Cap�tulo 10 do livro Telecommunications Breakdown.
clear
%-- pulrecsig.m - Criando a forma de pulso do sinal recebido.
N = 10000; M = 20 ; Ts = 0.0001;   % N� de s�mbolo, Fator de sobreamostragem
time = Ts*(N*M-1); t = 0:Ts:time;  % Tempo de amostragem e vetor de tempo
m=rand(N,1);                       % Cria sequ�ncia de s�mbolos aleat�rios
mup = zeros(1,N*M);                % ??
mup(1:M:end) = m;                  % ??
ps = hamming(M);                   % Janelamento com comprimento M
s = filter(ps,1,mup);              % Filtragem da informa��o
f0 = 1000; phoff = -1.0;           % Freq. da Portadora e a fase
c = cos(2*pi*f0*t+phoff);          % Portadora
rsc = s.*c;                        % Sinal Modulado (small carrier)
rlc = (s+1).*c;                    % Sinal Modulado (large carrier)
%----------------Aplicando a FFT--------------------%
fftrlc = fft(rlc);                 % Espectro de rlc
fftrsc = fft(rsc);                 % Espectro de rsc
[m,imax] = max(abs(fftrsc(1:end/2))); % pico m�ximo
ssf = (0:length(t))/(Ts*length(t));   % vetor de frequ�ncia
freqL = ssf(imax);                    % pico de frequ�ncia
phaseL = angle(fftrlc(imax));         % pico de fase 
% plot (ssf,abs([0 fftrlc]));
% plot (ssf,abs([0 fftrsc]));
%[m,imax] = max(abs(fftrlc(1:end/2)));

% pllpreprocess.m - Pr� processamento feito pelo m�todo quadr�tico e BPF
r = rsc;                         % Portadora suprimida
q = r.^2;                        % N�o linearidade quadr�tica
f1 = 500;                        % Frequ�ncia central do FPB
ff = [0 0.38 0.39 0.41 0.42 1];  % Ordem do filtro (5 steps)
fa = [0 0 1 1 0 0];
h = remez(f1,ff,fa);             % BPF via filtro remez
rp = filter(h,1,q);              % Resultado do processamento

% Recupera��o de frequ�ncia e fase "desconhecida" usando FFT
fftrBPF = fft(rp);                       % Espectro de rBPF
[m,imax] = max(abs(fftrBPF(1:end/2)));   % Freq de pico m�ximo
ssf = (0:length(rp))/(Ts*length(rp));    % vetor de frequ�ncia.
freqS = ssf(imax);                       % pico de frequ�ncia.
phasep = angle(fftrBPF(imax));           % pico de fase
[IR,f] = freqz(h,1,length(rp),1/Ts);     % Resposta em freq. do filtro
[mi,im] = min(abs(f-freqS));             % freq. em que ocorre o pico
phaseBPF = angle(IR(im));                % < que BPF na freq. de pico
phaseS = mod(phasep - phaseBPF,pi);      % estimador de fase

% pllsd.m - Rastramento de fase usando o SD
Ts = 1/10000;                   % Intervalo de tempo
time = 10;                      
t = 0:Ts:time-Ts;               % vetor de tempo
f0 = 1000;                       % freq. da portadora
phoff = -0.8;                   % fase da onda
%rp = cos(4*pi*f0*t + 2*phoff);  % sinal recebido
rp = cos(4*pi*f0*t + 2*phoff) + randn(1,length(t));
mu = 0.01;                     % tamanho do step
theta = zeros(1,length(t));     % m�todo de estimativa
fc = 1000;
for i=0:0.2:10
theta(1) = i;
f1 = 25;                        % Coeficiente
h = ones(1,f1)/f1;              % M�dia dos coeficientes
z = zeros(1,f1);                % buffer da m�dia
%for i = 0:-0.2:-4;
%theta(1) = i;
for k=1:length(t)-1;            % Algoritmo executado
    filtin(k) = (rp(k) - cos(4*pi*fc*t(k) + 2*theta(k)))*sin(4*pi*fc*t(k)+2*theta(k));
    z = [z(2:f1), filtin(k)];    % cont�m entradas anteriores
    theta(k+1) = theta(k) - mu*fliplr(h)*z'; % update = z convoluindo com h
end
plot(theta)
xlabel('time');
ylabel('phase offset');
title('Squared Difference Loop')
hold on
end
%for i = 0:-0.2:-4;
%theta(1) = i;

% pllconverge - Simulate Phase Locked Loop (PLL)
Ts = 1/10000;
time = 10;
t = Ts:Ts:time;
f0 = 1000;
phoff = -0.8;
%rp = cos(4*pi*f0*t + 2*phoff);
rp = cos(4*pi*f0*t + 2*phoff) + randn(1,length(t));
fl = 25;
ff = [0 0.01 0.02 1];
fa = [1 1 0 0];
h = remez(fl,ff,fa);        % LPF design
mu = 0.01;
fc = 1000;
theta = zeros(1,length(t)); % m�todo de estimativa
for i=0:0.2:10
theta(1) = i;
z = zeros(1,fl+1);          % z cont�m entradas anteriores
for k = 1:length(t)-1
    z = [z(2:fl+1), rp(k)*sin(4*pi*fc*t(k)+2*theta(k))];
    update = fliplr(h)*z';  % Nova sa�da do LPF
    theta(k+1) = theta(k) - mu*update; % update do algoritmo
end
plot(theta)
hold on
xlabel('time');
ylabel('phase offset');
title('PLL')
end

% costasloop.m - Simula��o do Costas Loop usando os par�metros
% de entrada do pulrecsig (gerado no in�cio do c�digo)
t = Ts:Ts:2*time;
%rp = cos(4*pi*f0*t + 2*phoff);
rp = cos(4*pi*f0*t + 2*phoff) + randn(1,length(t));
phoff = -0.8;
c = cos(2*pi*f0*t+phoff);
rsc = s.*c;
r = rsc;
fl = 25;
ff = [0 0.01 0.02 1];
fa = [1 1 0 0];
h = remez(fl,ff,fa);
mu = 0.01;
fc = 1000;
theta = zeros(1, length(t));
for i=0:0.2:10
theta(1) = i;
zs = zeros(1,fl+1);
zc = zeros(1,fl+1);
for k = 1:length(t)-1
    zs = [zs(2:fl+1), 2*r(k)*sin(2*pi*fc*t(k)+theta(k))];
    zc = [zc(2:fl+1), 2*r(k)*cos(2*pi*fc*t(k)+theta(k))];
    lpfs = fliplr(h)*zs';
    lpfc = fliplr(h)*zc';
    theta(k+1) = theta(k) - mu*lpfs*lpfc;
end
plot(theta)
hold on
xlabel('time');
ylabel('phase offset');
title('Costas Loop');
end

