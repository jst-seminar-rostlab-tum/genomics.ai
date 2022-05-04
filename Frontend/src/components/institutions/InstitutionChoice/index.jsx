import { Button } from '@mui/material';
import React, { useEffect, useState } from 'react';
import InstitutionService from 'shared/services/Institution.service';
import InstitutionList from '../InstitutionList';

function InstitutionChoice({ onChoiceChange }) {
  const [institutions, setInstitutions] = useState([]);
  const [chosen, setChosen] = useState(institutions[0]);

  useEffect(() => {
    InstitutionService.getMyAdminInstitutions()
      .then((newInstitutions) => {
        setInstitutions(newInstitutions);
        if (newInstitutions.length > 0) {
          setChosen(newInstitutions[0].id);
          onChoiceChange(newInstitutions[0].id);
        }
      })
      .catch(console.error);
  }, []);

  if (institutions.length === 0) {
    return <b>You are not an admin of any institution.</b>;
  }

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
          sx={{ color: institution.id === chosen ? 'green' : null }}
        >
          {institution.id === chosen ? 'Chosen' : 'Choose'}
        </Button>
      )}
    />
  );
}

export default InstitutionChoice;
