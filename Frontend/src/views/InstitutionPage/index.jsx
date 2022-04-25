import React, { useState, useEffect } from 'react';
import { useParams } from 'react-router-dom';
import TextField from '@mui/material/TextField';
import { getInstitution } from 'shared/services/mock/institutions';
import HeaderView from 'components/general/HeaderView';
import InstitutionMemberList from 'components/institutions/InstitutionMemberList';
import getUser from 'shared/services/mock/user';
import styles from './institutionPage.module.css';

function InstitutionPage() {
  let { id } = useParams();
  id = parseInt(id);
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

  useEffect(() => {
    getInstitution(id)
      .then((newInst) => { setTeam(newInst); updateIsAdmin(); })
      .catch((ignored) => { console.error(ignored); });
  }, [setInstitution, isAdmin]);

  useEffect(() => {
    getInstitution(id)
      .then((newInst) => { setInstitution(newInst); })
      .catch((ignored) => { console.error(ignored); });
  }, [setInstitution, isAdmin]);

  return (
    <>
    <div style={{ backgroundImage: `url(${institution.backgroundPictureURL})`, resizeMode: "center", flex: 1, height:400, width: undefined}}>
      <div class ={styles.container} ><img style={{width: 200, height: 200, borderRadius: 200 / 2}} src={institution.profilePictureURL} /></div>
    </div>
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
          fullwidth
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
