import React, { useState, useEffect } from 'react';
import {
  useParams
} from 'react-router-dom';
import HeaderView from 'components/general/HeaderView';
import { getTeam } from 'shared/services/mock/teams';

function TeamPage({ sidebarShown }) {
  const { id } = useParams();
  const [team, setTeam] = useState({});
  useEffect(() => {
    getTeam(id)
      .then((newTeam) => setTeam(newTeam))
      .catch((ignored) => { console.error(ignored); });
  }, [setTeam]);


  return (
    <HeaderView
      sidebarShown={sidebarShown}
      title={team.name}
    >
      {JSON.stringify(team)}
    </HeaderView>
  );
}

export default TeamPage;
