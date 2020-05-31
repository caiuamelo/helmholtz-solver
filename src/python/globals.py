USE_CPP_IMPLEMENTATION = False


def set_cpp(b):
    global USE_CPP_IMPLEMENTATION
    USE_CPP_IMPLEMENTATION = b


def choose_impl(f_a, f_b):
    global USE_CPP_IMPLEMENTATION
    if USE_CPP_IMPLEMENTATION:
        return f_a
    else:
        return f_b
