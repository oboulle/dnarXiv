
import hashlib

"""
use SHA256 hashing on sequences
"""


def base10_to_lower_base_convert(i: int, base: int) -> str:
    """
    return the number i written in a lower base
    """
    result = ""
    while i > 0:
            result = str(i % base) + result
            i = i // base
    return "".join(result)


def base16_to_baseX(base16_hash: str, baseX: int) -> str:
    """
    convert a string representing a number in base 16 to a string representing a number in any base X lower than 10
    """
    base10_hash = int(base16_hash, base=16)
    baseX_hash = base10_to_lower_base_convert(base10_hash, baseX)
    return baseX_hash


def hash_string_to_formated_base4(input_string: str, hash_size: int) -> str:
    """
    apply SHA256 hashing to a string to get a hash in base4 and of a specified size
    a basic hash is 64 long (base 16); so 128 long when reduced to base 4
    we count the needed number of basic hash to get to the wanted hash size / 2, because the size of the hash will be doubled in base 4
    then the input string is hashed with a "A" +k added at the end ( k is an incremented hash_number)
    the multiples hash are then concatenated to have a total hash of the wanted size /2 (rounded up to 64)
    the total hash is then reduced in base 4, and brought to the exact wanted size
    """
    hash_number = (hash_size/2) // 64 +1 
    total_hash = ""
    for k in range(int(hash_number)):
        # hash the string with a "A" and a number added, the A is to avoid getting the same hash for different strings
        # ex : hashing input = "a" with k = 10 equals hashing "a10"; but hashing input = "a1" with k = 0 is also equals hashing "a10"
        # here we will get the hashing of "aA10" and "a1A0", so completely different results
        # even if the total hashing cannot be the same, we try to avoid having some k_hash parts (of size 128) that are exactly the same between two different total hashes
        k_hash = sha256_hash(input_string + "A" + str(k))
        total_hash += k_hash
    base4_total_hash = base16_to_baseX(total_hash, 4)
    
    return base4_total_hash[:hash_size] # remove excess size


def hash_string_to_formated_base2(input_string: str, hash_size: int) -> str:
    """
    apply SHA256 hashing to a string to get a hash in base2 and of a specified size
    a basic hash is 64 long (base 16); so 256 long when reduced to base 2
    we count the needed number of basic hash to get to the wanted hash size / 2, because the size of the hash will be doubled in base 4
    then the input string is hashed with a "A" +k added at the end ( k is an incremented hash_number)
    the multiples hash are then concatenated to have a total hash of the wanted size /2 (rounded up to 64)
    the total hash is then reduced in base 4, and brought to the exact wanted size
    """
    hash_number = (hash_size/4) // 64 +1 
    total_hash = ""
    for k in range(int(hash_number)):
        # hash the string with a "A" and a number added, the A is to avoid getting the same hash for different strings
        # ex : hashing input = "a" with k = 10 equals hashing "a10"; but hashing input = "a1" with k = 0 is also equals hashing "a10"
        # here we will get the hashing of "aA10" and "a1A0", so completely different results
        # even if the total hashing cannot be the same, we try to avoid having some k_hash parts (of size 128) that are exactly the same between two different total hashes
        k_hash = sha256_hash(input_string + "A" + str(k))
        total_hash += k_hash
    base4_total_hash = base16_to_baseX(total_hash, 2)
    
    return base4_total_hash[:hash_size] # remove excess size


def sha256_hash(input_string: str) -> str:
    """
    SHA256 hashing, result is in base 16
    """
    h = hashlib.new('sha256')
    h.update(input_string.encode('utf-8'))

    return h.hexdigest()


# =================== main ======================= #
if __name__ == '__main__':
    for i in range(10):
        h = hash_string_to_formated_base4(str(i), 100)
        print(h, len(h))

