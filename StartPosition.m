function [ pos] =StartPosition( str, primer)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
k=strfind(str,primer);
pos=k+length(primer);

end

