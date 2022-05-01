import React, { useState, useEffect } from 'react';
import { useParams } from 'react-router-dom';
import TextField from '@mui/material/TextField';
import { getInstitution } from 'shared/services/mock/institutions';
import InstitutionMemberList from 'components/institutions/InstitutionMemberList';
import getUser from 'shared/services/mock/user';
import styles from './institutionPage.module.css';
import Avatar from '@mui/material/Avatar';

import queryMyTeams from 'shared/services/mock/teams';
import InstitutionTeamCard from 'components/institutions/InstitutionTeamCard';

function InstitutionPage() {
  let { id } = useParams();
  id = parseInt(id, 10);
  const [institution, setInstitution] = useState({});
  const [isAdmin, setIsAdmin] = useState(false);
  const [user, setUser] = useState({});

  const [teams, setTeams] = useState([]);
  useEffect(() => {
    queryMyTeams()
      .then((newTeams) => setTeams(newTeams))
      .catch((ignored) => { console.log(ignored); });
  }, [setTeams]);

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
      .then((newInst) => { setInstitution(newInst); updateIsAdmin(); })
      .catch((ignored) => { console.error(ignored); });
  }, [setInstitution, isAdmin]);

  useEffect(() => {
    getInstitution(id)
      .then((newInst) => { setInstitution(newInst); })
      .catch((ignored) => { console.error(ignored); });
  }, [setInstitution, isAdmin]);

  function onLeft(team) {
    setTeams(teams.filter((i) => i.id !== team.id));
  }

  return (
    <>
      <div className={styles.background} style={{ backgroundImage: `url(${institution.backgroundPictureURL})`, resizeMode: 'stretch' }}>
        <div className={styles.institutionIcon}>
          <Avatar src={institution.profilePictureURL} sx={{ width: 200, height: 200 }} />
        </div>
        <h1 className={styles.imageText}><span>{institution.name}</span></h1>
        <h3 className={styles.imageText}><span>{institution.country}</span></h3>
        <p className={styles.imageText}>
          <span>
            {institution.memberIds?.length + institution.adminIds?.length}
            Members
          </span>
        </p>
      </div>
      <div className={styles.test}>
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
            style={{ width: '100%' }}
            onChange={handleDescriptionChange}
            variant="standard"
          />
        </section>
        <section>
          <h2>Teams</h2>
          <hr />
          <div className={styles.content}>
            {teams.length === 0 ? 'No teams.' : ''}
            {teams.map((team) => (
              <div key={team.id}>
                <InstitutionTeamCard
                  team={team}
                  onLeft={(t) => onLeft(t)}
                />
                <div className={styles.cardSpacing} />
              </div>
            ))}
          </div>
        </section>
        <section>
          <h2>Members</h2>
          <hr />
          <InstitutionMemberList
            institution={institution}
            // eslint-disable-next-line no-shadow
            onRemoved={(institution, removedMember) => {
              setInstitution({
                ...institution,
                adminIds: institution.adminIds.filter((mId) => mId !== removedMember.id),
                memberIds: institution.memberIds.filter((mId) => mId !== removedMember.id),
              });
            }}
          />
        </section>
      </div>
    </>
  );
}

export default InstitutionPage;
