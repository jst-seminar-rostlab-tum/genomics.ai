import React, { useState, useEffect } from 'react';
import TextField from '@mui/material/TextField';
import getInstitutions from 'shared/services/mock/institutions';
import HeaderView from 'components/general/HeaderView';
import InstitutionMemberList from 'components/institutions/InstitutionMemberList';
import getUser from 'shared/services/mock/user';
//import styles from './institutionPage.module.css';

function InstitutionPage() {
  const [institution, setInstitution] = useState({});
  const [isAdmin, setIsAdmin] = useState(false);
  const [user, setUser] = useState({});

  function updateIsAdmin() {
    setIsAdmin((institution.adminIds || []).includes(user.id));
  }

  const handleDescriptionChange = (event) => {
    setInstitution({
      ...institution,
      description: event.target.value,
    });
  };

  useEffect(() => {
    getUser()
      .then((newUser) => { setUser(newUser); updateIsAdmin(); });
  }, [setUser, isAdmin]);

  return (
    <>
    <HeaderView>
    <section>
        <h2>Description</h2>
        <hr />
        <TextField
          id="description"
          multiline
          minRows={3}
          maxRows={5}
          value={institution.description}
          InputProps={{
            readOnly: !isAdmin,
          }}
          sx={{width: '100%'}}
          onChange={handleDescriptionChange}
          variant="standard"
        />
      </section>
      <section>
        <h2>Teams</h2>
        <hr />
      </section>
      <section>
        <h2>Members</h2>
        <hr />
        <InstitutionMemberList institution={institution} />
      </section>
      </HeaderView>
    </>
  );
}

export default InstitutionPage;
