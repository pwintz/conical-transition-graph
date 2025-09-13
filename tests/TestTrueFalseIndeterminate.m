classdef TestTrueFalseIndeterminate < pwintz.unittest.PWintzTestCase
% ! To run tests
% ` runtests TestTrueFalseIndeterminate

methods (Test)
  function test_canConvertToLogical(testClass)
    testClass.assertTrue(logical(TrueFalseIndeterminate.true));
    testClass.assertFalse(logical(TrueFalseIndeterminate.false));
    testClass.assertFalse(logical(TrueFalseIndeterminate.indeterminate));
  end % End of function.

  % ╭────────────────────────────────────╮
  % │             Truthiness             │
  % ╰────────────────────────────────────╯
  function test_truthiness_trueIsTruthy(testCase)
    % testCase.assertTru
    % testCase.assertTruthy(TrueFalseIndeterminate.true);
    if TrueFalseIndeterminate.true
      % OK
    else
      testCase.fatalAssertFail();
    end
  end % End of function.

  function test_truthiness_falseIsFalsy(testCase)
    % testCase.assertTru
    % testCase.assertTruthy(TrueFalseIndeterminate.true);
    if TrueFalseIndeterminate.false
      testCase.fatalAssertFail();
    end
  end % End of function.

  function test_truthiness_indeterminateIsFalsey(testCase)
    if TrueFalseIndeterminate.indeterminate
      testCase.fatalAssertFail();
    end
  end % End of function.

  % ╭────────────────────────────────────────╮
  % │             "And" operator             │
  % ╰────────────────────────────────────────╯
  function test_and_trueAndTrueIsTrue(testCase)
    left  = TrueFalseIndeterminate.true;
    right = TrueFalseIndeterminate.true;
    % ⋘──────── Vertify ────────⋙
    testCase.assertEqual(left  & right, TrueFalseIndeterminate.true);
    % Short circuit (uses conversion to logical)
    testCase.assertEqual(left && right, true);
  end % End of function.

  function test_and_trueAndFalseIsFalse(testCase)
    left  = TrueFalseIndeterminate.true;
    right = TrueFalseIndeterminate.false;
    % ⋘──────── Vertify ────────⋙
    testCase.assertEqual(left  & right, TrueFalseIndeterminate.false);
    % Short circuit (uses conversion to logical)
    testCase.assertEqual(left && right, false);
  end % End of function.
  
  function test_and_trueAndIndeterminateIsIndeterminate(testCase)
    left  = TrueFalseIndeterminate.true;
    right = TrueFalseIndeterminate.indeterminate;
    % ⋘──────── Vertify ────────⋙
    testCase.assertEqual(left  & right, TrueFalseIndeterminate.indeterminate);
    % Short circuit (uses conversion to logical)
    testCase.assertEqual(left && right, false);
  end % End of function.
  
  function test_and_falseAndIndeterminateIsFalse(testCase)
    left  = TrueFalseIndeterminate.false;
    right = TrueFalseIndeterminate.indeterminate;
    % ⋘──────── Vertify ────────⋙
    testCase.assertEqual(left  & right, TrueFalseIndeterminate.false);
    % Short circuit (uses conversion to logical)
    testCase.assertEqual(left && right, false);
  end % End of function.

  % ╭───────────────────────────────────────╮
  % │             "Or" operator             │
  % ╰───────────────────────────────────────╯
  function test_or_trueOrTrueIsTrue(testCase)
    left  = TrueFalseIndeterminate.true;
    right = TrueFalseIndeterminate.true;
    % ⋘──────── Vertify ────────⋙
    testCase.assertEqual(left  | right, TrueFalseIndeterminate.true);
    % Short circuit (uses conversion to logical)
    testCase.assertEqual(left || right, true);
  end % End of function.

  function test_or_trueOrFalseIsTrue(testCase)
    left  = TrueFalseIndeterminate.true;
    right = TrueFalseIndeterminate.false;
    % ⋘──────── Vertify ────────⋙
    testCase.assertEqual(left  | right, TrueFalseIndeterminate.true);
    % Short circuit (uses conversion to logical)
    testCase.assertEqual(left || right, true);
  end % End of function.
  
  function test_or_trueOrIndeterminateIsTrue(testCase)
    left  = TrueFalseIndeterminate.true;
    right = TrueFalseIndeterminate.indeterminate;
    % ⋘──────── Vertify ────────⋙
    testCase.assertEqual(left  | right, TrueFalseIndeterminate.true);
    % Short circuit (uses conversion to logical)
    testCase.assertEqual(left || right, true);
  end % End of function.
  
  function test_or_falseOrIndeterminateIsIndeterminate(testCase)
    left  = TrueFalseIndeterminate.false;
    right = TrueFalseIndeterminate.indeterminate;
    % ⋘──────── Vertify ────────⋙
    testCase.assertEqual(left  | right, TrueFalseIndeterminate.indeterminate);
    % Short circuit (uses conversion to logical)
    testCase.assertEqual(left || right, false);
  end % End of function.

  % ╭───────────────────────────────────────╮
  % │             "Not" operator            │
  % ╰───────────────────────────────────────╯
  function test_notTrueIsFalse(testCase)
    testCase.assertEqual(~TrueFalseIndeterminate.true, TrueFalseIndeterminate.false);
  end % End of function.

  function test_notIndeterminateIsIndeterminate(testCase)
    testCase.assertEqual(~TrueFalseIndeterminate.indeterminate, TrueFalseIndeterminate.indeterminate);
  end % End of function.
  
  function test_notFalseIsTrue(testCase)
    testCase.assertEqual(~TrueFalseIndeterminate.false, TrueFalseIndeterminate.true);
  end % End of function.
end % End of test methods

end % End of test methods block
