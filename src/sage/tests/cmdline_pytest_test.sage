# Test that pytest functions can also be defined in sage files

def simple_calculation():
    k.<a> = GF(5^3)
    assert a^124 == 1
