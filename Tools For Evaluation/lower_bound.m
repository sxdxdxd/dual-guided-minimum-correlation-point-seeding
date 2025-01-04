function Result = lower_bound(array, value)
first = 1;
last = length(array);
while first < last
    mid = first + floor((last - first) / 2);
    if array(mid) < value
        first = mid + 1 ;
    else
        last = mid ;
    end
end
Result = first;
end