function [rat,serial,bank] = decode_identifier(identifier)
    % parses a recording identifier into its component bits
    % input can be char or string, outputs are strings
    identifier = char(identifier);
    n = length(identifier);
    rat=string(identifier(1:4));
    serial=string(identifier(5:end-1));
    bank=string(identifier(end));
end