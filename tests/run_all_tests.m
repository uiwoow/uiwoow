% Run all fem2d framework tests and report pass/fail.
startup_fem2d;

tests = {
    'test_truss_bar',           @test_truss_bar;
    'test_beam_cantilever',     @test_beam_cantilever;
    'test_T3_patch',            @test_T3_patch;
    'test_Q4_patch',            @test_Q4_patch;
    'test_mixed_pure_beam',     @test_mixed_pure_beam;
    'test_mixed_pure_Q4',       @test_mixed_pure_Q4;
    'test_mixed_beam_to_Q4',    @test_mixed_beam_to_Q4;
};

n_pass = 0; n_fail = 0;
fprintf('\n========== fem2d test suite ==========\n');
for i = 1:size(tests, 1)
    name = tests{i,1};
    fn   = tests{i,2};
    fprintf('\n--- %s ---\n', name);
    try
        ok = fn();
        if ok
            n_pass = n_pass + 1;
        else
            n_fail = n_fail + 1;
            fprintf('  *** FAILED ***\n');
        end
    catch ME
        n_fail = n_fail + 1;
        fprintf('  *** ERROR: %s ***\n', ME.message);
    end
end
fprintf('\n========== Summary ==========\n');
fprintf('  Passed: %d / %d\n', n_pass, size(tests,1));
fprintf('  Failed: %d / %d\n', n_fail, size(tests,1));
if n_fail == 0
    disp('  ALL TESTS PASSED');
else
    disp('  SOME TESTS FAILED');
end
