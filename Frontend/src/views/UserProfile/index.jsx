/* eslint-disable react/jsx-props-no-spreading */
/* eslint-disable react/jsx-one-expression-per-line */
import React, { useState, useEffect } from 'react';
import HeaderView from 'components/general/HeaderView';
import InstitutionList from 'components/institutions/InstitutionList';
import TabPanel from 'components/TabPanel';
import Tabs from '@mui/material/Tabs';
import Tab from '@mui/material/Tab';
import Grid from '@mui/material/Grid';
import Avatar from '@mui/material/Avatar';
import TeamCard from 'components/teams/TeamCard';
import UserService from 'shared/services/User.service'
import stringToColor from 'shared/utils/stringColor';
import {
  useParams
} from "react-router-dom";
import CircularProgress from '@mui/material/CircularProgress';
import Box from '@mui/material/Box';

function a11yProps(index) {
  return {
    id: `simple-tab-${index}`,
    'aria-controls': `simple-tabpanel-${index}`,
  };
}

const testInstitutions = [
  {
    id: 5,
    name: 'Helmholtz Institute',
    country: 'Germany',
    profilePictureURL: 'https://www.hzdr.de/db/Pic?pOid=55058',
    backgroundPictureURL: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAAXNSR0IArs4c6QAAAA1JREFUGFdjsC9c9x8ABK8CXkzrf1cAAAAASUVORK5CYII=',
    adminIds: [],
    memberIds: [],
  },
  {
    id: 6,
    name: 'Technische Universität München',
    country: 'Germany',
    profilePictureURL: 'https://scalings.eu/wp-content/uploads/2019/11/tum-logo.png',
    backgroundPictureURL: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAAXNSR0IArs4c6QAAAA1JREFUGFdjsC9c9x8ABK8CXkzrf1cAAAAASUVORK5CYII=',
    adminIds: [1],
    memberIds: [],
  },
  {
    id: 7,
    name: 'Rostlab',
    country: 'Germany',
    profilePictureURL: 'https://avatars.githubusercontent.com/u/4093405?s=200&v=4',
    backgroundPictureURL: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAAXNSR0IArs4c6QAAAA1JREFUGFdjsC9c9x8ABK8CXkzrf1cAAAAASUVORK5CYII=',
    adminIds: [1],
    memberIds: [],
  },
];

const mockTeams = [
  {
    id: 1,
    name: 'Biotechnology Team',
    description: 'Biotechnology Team bla bla',
    adminIds: [6, 7],
    invitedMemberIds: [],
    memberIds: [1, 2, 3, 4, 5],
    visibility: 'public',
    institutionId: 2,
  },
];

function UserProfile({ sidebarShown }) {
  const [value, setValue] = useState(0);
  let { id } = useParams();


  const handleChange = (event, newValue) => {
    setValue(newValue);
  };

  const [user, setUser] = useState({})
  const [isUserLoading, setIsUserLoading] = useState(false)
  const [userInstitutions, setUserInstitutions] = useState([])
  const [isInstitutionLoading, setIsInstitutionLoading] = useState(false)
  const [userTeams, setUserTeams] = useState([])
  const [isTeamLoading, setIsTeamLoading] = useState(false)


  useEffect(() => {
    setIsUserLoading(true)
    UserService.getUser(id).then((data) => setUser(data)).finally(setIsUserLoading(false))
  }, [])

  useEffect(() => {
    setIsInstitutionLoading(true)
    UserService.getUserInstitutions(id).then((data) => setUserInstitutions(data)).finally(setIsInstitutionLoading(false))
  }, [])

  useEffect(() => {
    setIsTeamLoading(true)
    UserService.getUserTeams(id).then((data) => setUserTeams(data)).finally(setIsTeamLoading(false))
  }, [])

  return (
    <>
      <HeaderView
        sidebarShown={sidebarShown}
        title="Profile"
      >
        {isUserLoading ? <Box sx={{ display: 'flex' }}>
          <CircularProgress />
        </Box> :
          <Grid container justifyContent="center" alignItems="center" spacing={3}>
            <Grid item>
              <Avatar
                src={user.avatarUrl}
                alt={`${user.firstName} ${user.lastName}`}
                sx={{
                  backgroundColor: stringToColor(`${user.firstName} ${user.lastName}`),
                  width: 160,
                  height: 160,
                }}
              >
                {(user.firstName || '?')[0]}
              </Avatar>
            </Grid>
            <Grid item>
              <h2>
                {user.firstName} {user.lastName}
              </h2>
              <span>{user.email}</span>
            </Grid>
          </Grid>}
        <Box sx={{ width: '100%' }}>
          <Box sx={{ borderBottom: 1, borderColor: 'divider' }}>
            <Tabs value={value} onChange={handleChange} aria-label="basic tabs example">
              <Tab label="Institutions" {...a11yProps(0)} />
              <Tab label="Teams" {...a11yProps(1)} />
            </Tabs>
          </Box>
          <TabPanel value={value} index={0}>
            <InstitutionList
              isLoading={isInstitutionLoading}
              institutions={userInstitutions}
            />
          </TabPanel>
          <TabPanel value={value} index={1}>
            {userTeams.length === 0 ? 'No teams.' : ''}
            {userTeams.map((team) => (
              <div key={team.id}>
                <TeamCard
                  team={team}
                  onLeft={() => { }}
                />
              </div>
            ))}

          </TabPanel>
        </Box>
      </HeaderView>
    </>
  );
}

export default UserProfile;
