function seed = makeDeterministicSeed(taskId, neuronName, receptiveFieldType, foldIdentity)
%MAKEDETERMINISTICSEED Stable uint32 seed from task/cell/RF/fold identity.

key = [seedPart(taskId), '|', seedPart(neuronName), '|', ...
    seedPart(receptiveFieldType), '|', seedPart(foldIdentity)];

hash = uint32(2166136261);
prime = 16777619;
codes = uint16(char(key));

for i = 1:numel(codes)
    hash = bitxor(hash, uint32(codes(i)));
    hash = uint32(mod(double(hash) * prime, 2^32));
end

seed = double(hash);
if seed == 0
    seed = 1;
end
end

function part = seedPart(value)
if isnumeric(value)
    part = num2str(value);
elseif isstring(value)
    part = char(value);
elseif ischar(value)
    part = value;
else
    part = char(string(value));
end
end
