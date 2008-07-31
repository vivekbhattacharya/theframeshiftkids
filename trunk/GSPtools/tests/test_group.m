function test_group()
    result1 = group([1 2 2 2 2 2 3 4 4 4 4 4 5 6 7 8 0 0 0 2]);
    expected1 = [1 2 3 4 5 6 7 8 0 2; 1 5 1 5 1 1 1 1 3 1];
    fail_if(result1 ~= expected1, 'First');

    result2 = group([1 2 3 3 5 6 7 8 8 4 3]);
    expected2 = [1 2 3 5 6 7 8 4 3; 1 1 2 1 1 1 2 1 1];
    fail_if(result2 ~= expected2, 'Second');

    fail_if([] ~= group([]), 'Third');
    fail_if([3; 99] ~= group([3]), 'Fourth');
    disp 'All tests passed.';
end

function fail_if(boolean, which_one)
    if any(boolean)
        error(sprintf('%s test failed', which_one));
    end
end
