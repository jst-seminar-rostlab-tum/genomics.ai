import React, { useState, useEffect } from 'react';
import {
  useParams,
} from 'react-router-dom';
import HeaderView from 'components/general/HeaderView';
import { getTeam } from 'shared/services/mock/teams';
import getUser from 'shared/services/mock/user';
import { getInstitution } from 'shared/services/mock/institutions';

function TeamPage({ sidebarShown }) {
  const { id } = useParams();
  const [team, setTeam] = useState({});
  const [user, setUser] = useState({});
  const [institution, setInstitution] = useState({});
  const [isAdmin, setIsAdmin] = useState(false);

  function updateIsAdmin() {
    console.log(team, user);
    setIsAdmin((team.adminIds || []).indexOf(user.id) > -1);
  }

  useEffect(() => {
    getUser()
      .then((newUser) => { setUser(newUser); updateIsAdmin(); });
  }, [setUser, isAdmin]);

  useEffect(() => {
    getTeam(id)
      .then((newTeam) => { setTeam(newTeam); updateIsAdmin(); })
      .catch((ignored) => { console.error(ignored); });
  }, [setTeam, isAdmin]);

  // Institution may be undefined
  useEffect(() => {
    getInstitution(team.institutionId)
      .then((newInstitution) => setInstitution(newInstitution));
  }, [setInstitution]);

  return (
    <HeaderView
      sidebarShown={sidebarShown}
      title={team.name}
      rightOfTitle={<h4>{institution.name}</h4>}
    >
      {JSON.stringify(team)}
      {JSON.stringify(institution)}
      <br />
      {isAdmin ? 'Admin' : ''}
    </HeaderView>
  );
}

export default TeamPage;
