package ParseGenbank;
use strict; use warnings;
use Parse::RecDescent;
use Smooth;

my $grammar = q {
range: /\d+/ '..' /\d+/
{ $::selfish->range($item[1], $item[3]) }

join: 'join(' expr(s /,/) ')'
{ $::selfish->join(@{$item[2]}) }

complement: 'complement(' expr ')'
{ $::selfish->complement($item[2]) }

expr: join
    | complement
    | range
};

sub new {
    my ($class, $seq) = @_;

    my $self = {};
    $self->{parser} = new Parse::RecDescent($grammar);
    $self->{seq} = $seq;

    bless($self, $class); return $self;
}

sub parse {
    my ($self, $str) = @_;
    $self->{complemented} = 0;
    no warnings 'once';
    local $::selfish = $self;
    return $self->{parser}->expr($str);
}

sub range {
    my ($self, $a, $b) = @_;
    substr($self->{seq}, $a - 1, $b - $a + 1);
}

sub join {shift; join '', @_;}

sub complement {
    my $self = shift;
    $self->{complemented} = 1;
    Smooth::reverse_complement shift;
}

if ($0 eq __FILE__) {
    use Data::Dumper;
    use 5.010;
    $Data::Dumper::Indent = 0;
    my $parser = new ParseGenbank('aucgaucg');
    say Dumper $parser->parse('complement(join(1..8,1..8,1..8))');
    say Dumper $parser->parse('complement(1..8)');
    say Dumper $parser->parse('join(complement(join(1..2,2..3)),5..6)');
    say Dumper $parser->parse('join(1..2,2..3)');
    say Dumper $parser->parse('1..2');
}

1;
