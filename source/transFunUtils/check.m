function check(message,condition)

  if condition
    str = ['Checking ' message '... [PASSED]'];
  else
    str = ['Checking ' message '... [FAILED]'];
  end

  disp(str);

end
