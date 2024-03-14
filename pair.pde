public class Key {
  public final int first;
  public final int second;
  public Key(final int _f, final int _s) {
    this.first = _f;
    this.second = _s;
  }
  @Override
  public boolean equals (final Object O) {
    if (!(O instanceof Key)) return false;
    if (((Key) O).first != first) return false;
    if (((Key) O).second != second) return false;
    return true;
  }
  @Override
  public int hashCode() {
    //return (first << 16) + second;
    return Objects.hash(first, second);
  }
  public Key opp() {
    return new Key(second, first);
  }
  public Key sorted() {
    if (first <= second) return this;
    return opp();
  }
}

public class SortedKey extends Key {
  public SortedKey(int _f, int _s) {
    super(_f <= _s ? _f : _s, _f <= _s ? _s: _f);
  }
}
