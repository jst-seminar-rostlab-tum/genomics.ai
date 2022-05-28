import Button from 'components/CustomButton';
import React, { useState } from 'react';
import InstitutionList from '../InstitutionList';
import { useInstitutions } from 'shared/context/institutionContext';

function InstitutionChoice({ onChoiceChange }) {
  const { institutions } = useInstitutions();
  const [chosen, setChosen] = useState(institutions[0]);

  return (
    <InstitutionList
      institutions={institutions}
      disableLinks
      replaceTrailing={(institution) => (
        <Button
          onClick={() => {
            setChosen(institution.id);
            onChoiceChange(institution.id);
          }}
          type={institution.id === chosen ? 'primary' : 'tertiary'}
        >
          {institution.id === chosen ? 'Chosen' : 'Choose'}
        </Button>
      )}
    />
  );
}

export default InstitutionChoice;
