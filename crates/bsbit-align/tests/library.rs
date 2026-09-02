#![allow(missing_docs)]

use bsbit_align::library::{
    LibraryProfile, PairConstraintError, PairConstraints, TemplateSpan, TemplateSpanBounds,
};

#[test]
fn paired_library_values_validate_without_search_policy() {
    let bounds = TemplateSpanBounds::new(TemplateSpan::new(7), TemplateSpan::new(31))
        .expect("ordered bounds");
    assert!(bounds.contains(TemplateSpan::new(7)));
    assert!(bounds.contains(TemplateSpan::new(31)));
    assert!(!bounds.contains(TemplateSpan::new(32)));

    let constraints = PairConstraints::new(LibraryProfile::Directional, bounds);
    assert_eq!(constraints.profile(), LibraryProfile::Directional);
    assert_eq!(constraints.span_bounds(), bounds);

    assert_eq!(
        TemplateSpanBounds::new(TemplateSpan::new(32), TemplateSpan::new(31)),
        Err(PairConstraintError::InvertedSpanBounds {
            minimum: TemplateSpan::new(32),
            maximum: TemplateSpan::new(31),
        })
    );
}
