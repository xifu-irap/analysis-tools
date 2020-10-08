
###############################################################################
def switch_bin2hexa_digit(digit):
    switcher={
    "0000": "0",
    "0001": "1",
    "0010": "2",  
    "0011": "3", 
    "0100": "4",   
    "0101": "5",   
    "0110": "6",   
    "0111": "7",   
    "1000": "8",   
    "1001": "9",   
    "1010": "A",   
    "1011": "B",   
    "1100": "C",   
    "1101": "D",   
    "1110": "E",   
    "1111": "F"  
    }
    return(switcher[digit])


###############################################################################
def twos_comp(val, bits):
    r"""
    This function computes the 2's complement of int value val
    """
    if (val & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
        val = val - (1 << bits)        # compute negative value
    return val                         # return positive value as is

###############################################################################
def twos_comp_bin2int(binary_string):
    r"""
    This function converts a two's complement binary value (string) to an integer
    """
    return(twos_comp(int(binary_string,2), len(binary_string)))

###############################################################################
def twos_comp_hex2int(hex_string):
    r"""
    This function converts a two's complement hex value (string) to an integer
    """
    return(twos_comp(int(hex_string,16), 32))

###############################################################################
def switch_bin2hexa(bin):
    r"""
    This function makes the conversion from binary ascii to hexadecimal ascci.
    """

    nbits=len(bin)
    offset=0
    hexa=''
    while offset<nbits:
        hexa=hexa+switch_bin2hexa_digit(bin[offset:offset+4])
        offset+=4
    return(hexa)

###############################################################################
def switch_natbin2dec(bin):
    r"""
    This function makes the conversion from natural binary ascii to decimal ascci.
    """

    nbits=len(bin)
    offset=0
    dec=''
    while offset<nbits:
        dec=dec+bin[offset]*2**offset
        offset+=1
    return(dec)

###############################################################################
def dec2cad(dec, n):
    r"""
    This function computes the 2s complement binary value of an integer.
    
    Parameters
    ----------
    dec : number
        The decimal input integer to be converted
    n : number
        The number of bits of the output string

    Output
    ------
    bin_str : string
        The binary value in a string format

    Raises
    ------
    ValueError
        When n is too small to do the conversion

    See Also
    --------
    num : string to number conversion

    Examples
    --------
    >>> dec2cad(7, 8)
    '00000111'

    >>> dec2cad(-7, 8)
    '11111001'

    >>> dec2cad(7, 2)
    ValueError: Requested size is too small for this value

    """
    bin_str = ''
    if dec >= 0:
        while n > 0:
            bin_str = str(dec % 2) + bin_str
            dec >>= 1
            n += -1
    else:
        dec = abs(dec + 1)
        while n > 0:
            rev_str = ('1', '0')
            bin_str = rev_str[dec % 2] + bin_str
            dec >>= 1
            n += -1
    if dec > 0:
        raise ValueError('Requested size is too small for this value')
    return bin_str

###############################################################################
def dec2natbin(dec, n):
    r"""
    This function computes the natural binary value of a positive integer.
    
    Parameters
    ----------
    dec : number
        The decimal input integer to be converted
    n : number
        The number of bits of the output string

    Output
    ------
    bin_str : string
        The binary value in a string format

    Raises
    ------
    ValueError
        When n is too small to do the conversion

    See Also
    --------
    num : string to number conversion

    """
    bin_str = ''
    if dec >= 0:
        while n > 0:
            bin_str = str(dec % 2) + bin_str
            dec >>= 1
            n += -1
    if dec > 0:
        raise ValueError('Requested size is too small for this value')
    return bin_str

###############################################################################

def signed_hexa_2_dec(hexstr,nbits):
    r"""
    This function computes the decimal value of a 2's complement signed hexa string.
    
    Parameters
    ----------
    hexstr : string
        The hexadecimal value to be converted
    nbits : number
        The number of bits of the data

    Output
    ------
    value : integer
        The converted value

    Examples
    --------
    >>> signed_hexa_2_dec('FFFE',16)
    -2
    >>> signed_hexa_2_dec('7FFF',16)
    32767
    >>> signed_hexa_2_dec('7F',8)
    127
    >>> signed_hexa_2_dec('FF',8)
    -1
    """
    value = int(hexstr,16)
    if value & (1 << (nbits-1)):
        value -= 1 << nbits
    return value

###############################################################################

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

    return(hex(signeddec & (2**nbits-1)))

###############################################################################
