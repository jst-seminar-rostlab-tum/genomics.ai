import React, { useState, useEffect } from 'react';
import {
  useParams,
} from 'react-router-dom';
import Chip from '@mui/material/Chip';
import Stack from '@mui/material/Stack';
import HeaderView from 'components/general/HeaderView';
import JobList from 'components/teams/detail/JobList';
import { getTeam } from 'shared/services/mock/teams';
import getUser from 'shared/services/mock/user';
import { getInstitution, queryIsAdminInstitutions } from 'shared/services/mock/institutions';

export default function TeamPage({ sidebarShown }) {
  const { id } = useParams();
  const [team, setTeam] = useState({});
  const [user, setUser] = useState({});
  const [institution, setInstitution] = useState({});
  const [isMember, setIsMember] = useState(false);
  const [isAdmin, setIsAdmin] = useState(false);
  const [adminInstitutions, setAdminInstitutions] = useState([]);

  function updateIsMember() {
    setIsMember((team.memberIds || []).indexOf(user.id) > -1);
  }

  function updateIsAdmin() {
    setIsAdmin((team.adminIds || []).indexOf(user.id) > -1);
  }

  useEffect(() => {
    getUser()
      .then((newUser) => { setUser(newUser); updateIsAdmin(); updateIsMember(); });
  }, [setUser, isAdmin, isMember]);

  useEffect(() => {
    getTeam(id)
      .then((newTeam) => { setTeam(newTeam); updateIsAdmin(); updateIsMember(); })
      .catch((ignored) => { console.error(ignored); });
  }, [setTeam, isAdmin, isMember]);

  // Institution may be undefined
  useEffect(() => {
    getInstitution(team.institutionId)
      .then((newInstitution) => setInstitution(newInstitution));
  }, [team, setInstitution]);

  useEffect(() => {
    queryIsAdminInstitutions(user.id)
      .then((newAdminInstitutions) => setAdminInstitutions(newAdminInstitutions));
  }, [user, setAdminInstitutions]);

  return (
    <HeaderView
      sidebarShown={sidebarShown}
      title={team.name}
      rightOfTitle={(
        <HeaderOptions
          team={team}
          isAdmin={isAdmin}
          institution={institution}
          isMember={isMember}
          availableInstitutions={adminInstitutions}
          setInstitution={setInstitution}
        />
      )}
    >
      <br />
      <section>
        <h2>Description</h2>
        <hr />
        TODO: implement
      </section>
      <section>
        <h2>GeneMapper</h2>
        <hr />
        <JobList teamId={id} forPart="geneMapper" />
      </section>

    </HeaderView>
  );
}

function HeaderOptions({
  team, isAdmin, institution,
}) {
  return (
    <Stack
      direction="row"
      alignItems="flex-end"
      spacing={2}
      sx={{ paddingLeft: '20px' }}
    >
      {(!isAdmin && institution) && <h4>{institution.name}</h4>}
      {!isAdmin && <Chip label={team.visibility} color="primary" />}

    </Stack>
  );
}
