from gam.converters.convert_three_letter_mutation import convert_three_letter_mutation


def test_convert_three_letter_mutation_standard_case() -> None:
    assert convert_three_letter_mutation("Ser450Leu") == "S450L"


def test_convert_three_letter_mutation_stop_case() -> None:
    assert convert_three_letter_mutation("Tyr103*") == "Y103!"
