import numpy as np

#------------------------------------------------------------
def to_signed(val, nbits):
    """This function converts a positive integer value into 
    a string corresponding to the signed hexa decimal value

    Args:
        val (array): values to be converted
        nbits (integer): number of bits of the value (position of the sign bit)

    Raises:
        ValueError: error raised in case of over range.
        
    Examples:
    >>> to_signed(7,4)
    '7'
    >>> to_signed(7,3)
    '-1'
    >>> to_signed(8,4)
    '-8'
    """    
    if max(val) >= 2**nbits:
        raise ValueError('There is one (or more) over-range on this {0:d}-bit data!'.format(nbits))
    
    # When sign bit is equal to one shifting to negative values 
    val[check_bit_is_one(val, nbits)] -= 2**nbits
    
    return(val)

#------------------------------------------------------------
def check_bit_is_one(data, bit_nb):
    return((data % (2**bit_nb)) // 2**(bit_nb-1) == 1)

#------------------------------------------------------------
def dec_to_signed_hexa(signeddec,nbits):
    r"""
    This function computes the hexadecimal string of a signed decimal value.
    
    Parameters
    ----------
    signeddec : decimal
        The decimal value to be converted
    nbits : number
        The number of bits of the data

    Output
    ------
    value : string
        The converted value

    Examples
    --------
    >>> dec_to_signed_hexa(2,8)                                                           
    '0x2'
    >>> dec_to_signed_hexa(2,16)                                                          
    '0x2'
    >>> dec_to_signed_hexa(-2,16)                                                         
    '0xfffe'
    """

    return(hex(signeddec & (2**nbits-1))[2:])

#------------------------------------------------------------
