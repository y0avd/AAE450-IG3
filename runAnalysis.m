
tic
launchLow = 2461771.5;
launchHigh = 2466519.5;
launch = linspace(launchLow,launchHigh,90);
TOF = linspace(6,18,90);
delV = zeros(1,length(TOF));
result = zeros(length(launch),length(TOF));

for k = 1:length(launch)
    for i = 1:length(TOF)
        delV(i) = calcDelV(launch(k),TOF(i));
    end
result(k,:) = delV(:);
end
mindV = 100;
figure(1)
for i = 1:length(launch)
    plot(TOF,result(i,:),'b')
    hold on
    for k = 1:length(TOF)
        if result(i,k) < mindV
            mindV = result(i,k);
            mindVTOF = TOF(k);
            mindVLaunch = launch(i);
        end
    end
end

weightedMin = zeros(length(launch),length(TOF));
min = 1000;
for i = 1:length(launch)
for k = 1:length(TOF)
    weightedMin(i,k) = 0.3*TOF(k)+0.7*result(i,k);
        if weightedMin(i,k) < min
            minTOFdV = result(i,k);
            minTOF = TOF(k);
            minTOFLaunch = launch(i);
            min = weightedMin(i,k);
        end
end
end

datedV = convertCharsToStrings(datestr(datetime(mindVLaunch,'convertfrom','juliandate')));
dateTOF = convertCharsToStrings(datestr(datetime(minTOFLaunch,'convertfrom','juliandate')));
dataLabeldV = sprintf("DeltaV: %0.3f km/s\nTOF: %0.3f years\nLaunch Date: %s",mindV,mindVTOF,datedV);
dataLabelTOF = sprintf("DeltaV: %0.3f km/s\nTOF: %0.3f years\nLaunch Date: %s",minTOFdV,minTOF,dateTOF);
hold on
plot(mindVTOF,mindV,'or','MarkerSize',8)
plot(minTOF,minTOFdV,'og','MarkerSize',8)
xlabel("TOF, years")
ylabel("Mission \DeltaV")
text(mindVTOF + 0.1, mindV - 3.5, dataLabeldV, 'BackgroundColor',[0.95 0.95 0.95])
text(minTOF + 0.1, minTOFdV - 4, dataLabelTOF, 'BackgroundColor',[0.95 0.95 0.95])
title('Earth Direct to Uranus Mission Frontier: Launch Dates 2028-01-01 to 2040-12-31')
grid on
toc