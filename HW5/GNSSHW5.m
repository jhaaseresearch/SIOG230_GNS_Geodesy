% Loading observation file via Jennifer's function:
[obstypes, t, gpssat, sitepos, hant] = read_rinexo('opmt2920.19o');
[eph,alpha,beta,dow] = read_rinexn('brdc2920.19n')
%Establishing Initial Variables:
c = 0.299792458 * 10^9;
initialX = 4202777.6220   
initialY = 171367.7610  
initialZ = 4778659.9840
X = [20,20,20,20]
%C1:
C1 = obstypes.C1;
%P1:
P1 = obstypes.P1;
%P2:
P2 = obstypes.P2;
tinsecs = (0:30:(30* 2880 -1)).';
tinsecs = tinsecs + (3600 * 24 * 6)
    usedtimes = []
    finalX = []
    finalY = []
    finalZ = []
    finallong = []
    finallat = []
    finalheight = []
    finalclockbiasrec = []
    CovarianceMatrix = []
    CovarianceMatrixlatlong = []
    VDOP = []
    HDOP = []
    PDOP = []
    TDOP = []
 for i = 1: length(tinsecs)
    t = tinsecs(i)
    row = C1(i, :)  
    A = []
    satxyz = []
    usedsats = []
    satclockbias = []
    usedobs = []
    geometricdist = []
    if sum(isnan(row)) <= 26 % Sum row nans > 26
    % Times for getsatpos:
    %usedtimes(end+1) = t(i)
        t_satpos = t - (C1(i,:)./c); %change row value
        usedtimes(end + 1) = t;
    %For all satellites get_satpos:
        counta = 0
        for k = 1:30
            if (isnan(row(k))) == 0
                counta = counta + 1
                [pos]  = get_satpos(t_satpos(1,k), gpssat(k), eph,0)%NOTE: I took a working get_satpos from Jennifer
                satxyz(1,counta) =pos(1)
                satxyz(2, counta)=pos(2)
                satxyz(3,counta) =pos(3)
                satclockbias(1,counta) = pos(4)
                usedsats(counta) = gpssat(k)
                usedobs(counta) = C1(i,k) %REPLACE WITH TIME INDEX
            end
        end
        littlel= []
        ax = []
        ay = []
        az = []
        ac = []
        for ind = 1: length(usedobs)
                geometricdist(ind) = sqrt((satxyz(1,ind) -initialX)^2 + (satxyz(2,ind)-initialY)^2 + (satxyz(3,ind)-initialZ).^2)
                littlel(ind) = usedobs(ind) - geometricdist(ind) + c*satclockbias(ind)
                ax(ind) = -(satxyz(1,ind) - initialX) / geometricdist(ind)
                ay(ind) = -(satxyz(2,ind)-initialY) / geometricdist(ind)
                az(ind) = -(satxyz(3,ind)-initialZ)/ geometricdist(ind)
        end
        for o = 1: length(ax)
            ac(o) = -c / 10^9
        end 
        A(:,1) = ax;
        A(:,2) = ay
        A(:,3) = az
        A(:,4) = ac

        while X(1) > 1 || X(2) > 1
            disp(A)
            C = A.' * A;
            X = (pinv(C)*A.')*(littlel.')
            guessX = initialX + X(1)
            initialX = guessX
            guessY = initialY + X(2)
            initialY = guessY
            guessZ = X(3)+initialZ
            initialZ = guessZ
            guessrecclockbias = X(4)
            littlel = []
            geometricdist= []
            for ind = 1: length(usedsats)
                geometricdist(ind) = sqrt((satxyz(1,ind) -guessX)^2 + (satxyz(2,ind)-guessY)^2 + (satxyz(3,ind)-guessZ)^2)
                littlel(ind) = usedobs(ind) - geometricdist(ind) + c*satclockbias(ind)
                ax(ind) = -(satxyz(1,ind) - guessX) / geometricdist(ind)
                ay(ind) = -(satxyz(2,ind)-guessY) / geometricdist(ind)
                az(ind) = -(satxyz(3,ind)-guessZ)/ geometricdist(ind)
            end
            A(:,1) = ax
            A(:,2) = ay
            A(:,3) = az
            A(:,4) = ac
        end
    disp(guessX)
    finalX(end + 1) = guessX
    finalY(end+1) = guessY
    finalZ(end+1) = guessZ
    finalclockbiasrec(end+1) = guessrecclockbias
    input = [t, finalX,finalY,finalZ]
    [time, long, lat, height] = xyz2wgs(input)
    finallong(end + 1) = long
    finallat(end + 1) = lat
    finalheight(end + 1) = height
    HELP = pinv(A.' *A)
    CovarianceMatrix(end + 1,:,:) = pinv(A.' *A)
    COVMAT= reshape(CovarianceMatrix(end,1:4,1:4),[4,4])
    lat = deg2rad(lat)
    long = deg2rad(long)
    R = [-sin(lat)*cos(long) -sin(lat)*sin(long) cos(lat)
    -sin(long) cos(long) 0
    cos(lat)*cos(long) cos(lat)*cos(long) sin(lat)]
    CovarianceMatrixlatlong(end +1,:,:) = R * COVMAT(1:3,1:3) * R.'
    THISMATRIX = reshape(CovarianceMatrixlatlong(end,1:3,1:3), [3,3])
    VDOP(end +1) = sqrt(THISMATRIX(3,3))
    HDOP(end +1) =sqrt(THISMATRIX(1,1)+ THISMATRIX(2,2))
    PDOP(end+1) = sqrt(trace(THISMATRIX))
    TDOP(end +1) = sqrt(COVMAT(4,4))
end

end
%%
figure(1)
subplot(3,1,1)
plot(usedtimes, finalX, '-', 'Color',[0 0 0], 'LineWidth',3)
xlabel('time in seconds')
ylabel('X position')
subplot(3,1,2)
plot(usedtimes, finalY, '-', 'Color',[0 0 0], 'LineWidth',3)
xlabel('time in seconds')
ylabel('Y position')
subplot(3,1,3)
plot(usedtimes, finalZ, '-', 'Color',[0 0 0], 'LineWidth',3)
xlabel('time in seconds')
ylabel('Z position')


