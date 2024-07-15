from TD2.translator import Translator

def test_translate_short():
    t = Translator()
    assert t.translate("ATG") == ("M", [0])

