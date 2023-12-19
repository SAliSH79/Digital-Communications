
%% Exp 4
% Teacher : MS_Jafari
% Author: [SeyedAli] - [SeyedHosseini]
% E-mail: [alishosseini79@aut.ac.ir]
%Student-Number : [9723042]
% University: Amirkabir University of Technology
%% Transmitter
    clc;
    close all;
    clear;
    %% Initialization Data
    clear
    clc
    s  = importdata('s_k.jpg');
    figure(1)
    imshow(s)
    grid on;
    title(" Image for test")

    imbins = zeros(size(s,1)*size(s,2)*3,8);
    for i = 1:size(s,1)
        for j = 1:size(s,2)
            for k = 1:size(s,3)
                imbins(k+3*((i-1)*size(s,2)+(j-1)),1:8)=de2bi(s(i,j,k),8,'left-msb');
            end
        end
    end
    Stream = reshape(imbins',1,3*8*size(s,1)*size(s,2));

    %% Channel Coding

    r = 5 ;
    n = 2^r - 1; %Codeword Length
    k = n - r; %Massage Length

    Encod_S = encode(Stream,n,k,'hamming/binary');

    L_Encd = length(Encod_S);
    L_UnEncd = length(Stream);

    %% Interleaving
    clc;
    N_Col = 31 ;
    N_Row = ceil(L_Encd/N_Col);
    InterleavedStream = reshape(reshape(Encod_S,N_Row,N_Col)' , 1 , N_Row*N_Col);

    figure(2)
    plot(abs(xcorr(InterleavedStream,Encod_S)))
    grid on;
    title(" Correlation of Encoded and Interleaved Signal")

    %% Bit to Integer
    clc;
    h = [+1 ,+1 ,+1 ,+1 ,+1 ,0 ,0, +1 ,+1 ,0 ,+1 ,0, +1];

    Header = repmat(h,1,5);
    clc;
    figure(3)
    plot(abs(xcorr(Header)),"k")
    grid on;
    ylabel("Mag")
    xlabel("time")
    title("AutoCorr Header")
    % axis([0 260 -5 140])

    Trans_Stream = [Header InterleavedStream];
    T_L = length(Trans_Stream);
    %% Initialization
    Rs = 1e4;

    SPS = 8;
    alpha = 0.2;
    Span = 10;
    %% part 1
    M = 16;
    N = 1e3; 
    D_T= Trans_Stream ;
    Data_App= [D_T,zeros(1,11091)];

    %% Binary to Decimal
    clc;
    Data_Shaped = reshape(Data_App,4,120000)';
    IntegerData = bi2de(Data_Shaped,'left-msb')';

    %% Modulation
    clc;
    % D = D - 1 ;
    D_QAM = qammod(IntegerData,M);

    H = [0,1,1,1,0,1,0,0,1,1,1,0,1,1,0,0,1,1,0,0,0,...
        0,1,1,1,0,0,1,1,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,1....
        ,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,1,0,1,0,0,1,1,1,....
        0,0,0,0,1,1,0,1,1,0,1,1,1,1,1,1,1,1,0,0,1,0,1....
        ,0,0,0,0,0,0,0,1,0,0,0,1,1,0,1,1,1,0,1,1,1,0,1....
        ,0,1,0,1,0,0,1,1,1,0,1,0,1,0];

    hh = numel(H);
    %% Plotting
    clc;
    figure(4)
    plot(abs(xcorr((-1).^H)),"r")
    grid on;
    ylabel("Mag")
    xlabel("time")
    title("AutoCorr Header of Phy")
    axis([0 260 -5 140])

    %% Pilot
    Ds = D_QAM;
    clc;
    P = (1 + 1i)/(sqrt(2));
    pp=20;
    pilot = repmat(P,1,pp);
    header = (-1).^H;

    %% Frame
    clc;
    % frame = [header Data pilot Data pilot Data];
    % % clc;
    % % f1 = [pilot Ds];
    % % f11 = repmat(f1,1,2);
    % % f22 = [header Ds];
    % % frame = [f22 f11]; %frame 
    %  L = length(frame);
    frame = zeros(40,3168);
    n = 1e3;
    for i = 1 : 40
        D11 = Ds((i-1)*n + 1 : i*n);
        D22 = Ds(i*n + 1 : (i+1)*n);
        D33 = Ds((i+1)*n + 1 : (i+2)*n);

        D11 = D11/3.1225;
        D22 = D22/3.1225;
        D33 = D33/3.1225;

        frame(i,:) = [header D11 pilot D22 pilot D33];
    end
%% Scatter Plot
    clc;
    Frm = reshape(frame',1,40*3168)';
    L = length(Frm);
    figure()
    scatterplot(Frm)
    grid on;
    ylabel("imag")
    xlabel("real")
    title("Constl 16 QAM with Header Pilot Transmitter")
    % axis([0 260 -5 140])


%% Receiver
    clc;

    %% Integer to Bit
    clc;
    x = out.yout;
    %Rec_Bits_Shaped = de2bi(IntegerData...
        %,'left-msb');
    Rec_bit = de2bi(x ,'left-msb') ;
    isequal(Rec_bit,Data_Shaped)
    %Rec_Data_App = reshape(Rec_Bits_Shaped'...
      %  ,1,4*120000);
    Rec_Data_App = reshape(Rec_bit'...
        ,1,4*662966);
    isequal(Rec_Data_App,Data_App)

    %% Filtering 
    Filtered_RecD = FilterD((-1).^Rec_Data_App,...
        (-1).^Header);
    figure()
    plot(Filtered_RecD)
    %% Detecting --- May find the error 
    clc;
    [Max,I] = maxk(Filtered_RecD,5); %find max and index
    LengthI = length(InterleavedStream);
    Rec_I_S = Rec_Data_App(I(4) + 1 : I(4) + LengthI);
    
    isequal(Rec_I_S,InterleavedStream)
    
    %% DeInterLeaving
    clc;
    Rec_Enc_DI = reshape(reshape(Rec_I_S,N_Col,N_Row)' , 1 , N_Row*N_Col);
    isequal(Rec_Enc_DI,Encod_S)
    
    %% Decoder
    clc;
    r = 5 ;
    n = 2^r - 1; %Codeword Length
    k = n - r; %Massage Length
    
    Decod_S = decode(Rec_Enc_DI , n, k, 'hamming/binary');
    Decod_S = Decod_S(1 : end - 8);
    isequal(Stream,Decod_S)
    
%% Showing Image
    %% Image Reconstruction
    StackedBins   = reshape(Decod_S,8,3*128^2)';
    RxImage = zeros(128,128,3);
    for i = 1:size(RxImage,1)
        for j = 1:size(RxImage,2)
            for k = 1:size(RxImage,3)
                RxImage(i,j,k) = bi2de(StackedBins(k+3*((i-1)*size(RxImage,2)+(j-1)),1:8),'left-msb');
            end
        end
    end
    RxImage = uint8(RxImage);
    %% Plotting
    figure()
    imshow(RxImage)
    grid on;
    title("Received Image")
   %% Plotting xcor
   figure()
   plot(abs(xcorr(Encod_S,Rec_Enc_DI)))