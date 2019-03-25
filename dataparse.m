%Written by Vince Chiu
%Takes the data that is spit out by cosmed and turns it into vO2 and vCO2

function [vo2,vco2,ready] = dataparse(Data)
%each row is a new instance of '<'
%each column contain the decimal characters up to the next '>'
% < = 60
% > = 62
maxnode = 40;
maxmessagelen = 60;
message = zeros(maxnode,maxmessagelen);
node = 0;
charidx = 1;
vo2 = 0;
vco2 = 0;
ready = 0;
%parse message
for msgidx = 1:length(Data)
    if charidx == maxmessagelen+1
        charidx = maxmessagelen; %prevent character overflow
    end
    if Data(msgidx) == 60 %<
        node = node + 1;
        charidx = 1;
    end
    message(node,charidx) = Data(msgidx);
    charidx = charidx + 1;
end
%find vo2 and vco2 nodes
for nodeidx = 1:maxnode
    temp = message(nodeidx,:);
    if all(temp(1:5) == [60 86 79 50 62]) %<VO2>
        vo2 = extractvalue(temp(6:end));
    elseif all(temp(1:6) == [60 86 67 79 50 62]) %<VCO2>
        vco2 = extractvalue(temp(7:end));
    elseif all(temp(1:8) == [60 82 101 115 117 108 116 62]) %<Result>
        if all(temp(9:11) == [65 67 75]) %ACK
            ready = 1;
        elseif all(temp(9:11) == [78 65 75]) %NAK
            ready = -1;
        end
    end
end

function value = extractvalue(data)
periodidx = find(data==46,1);
data = data - 48; %makes decimal numbers into ascii numbers
if periodidx == 4 %three digit vo2 or vco2
    value = data(1)*100 + data(2)*10 + data(3);
elseif periodidx == 5 %four digit vo2 or vco2
    value = data(1)*1000 + data(2)*100 + data(3)*10 + data(4);
else
    value = -1; %error
end