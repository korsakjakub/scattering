package main

import (
	"github.com/spf13/viper"
	"log"
)

type Config struct {
	Ratio           string `mapstructure:"ratio"`
	InitialPosition string `mapstructure:"initialposition"`
	InitialPhi      string `mapstructure:"initialphi"`
	CollisionEnergy string `mapstructure:"collisionenergy"`
	IndexRange      string `mapstructure:"indexrange"`
}

var vp *viper.Viper

func LoadConfig(additionalPath []string, args ...string) Config {
	vp = viper.New()
	var config Config
	if len(args) > 0 {
		vp.SetConfigName(args[0])
		vp.SetConfigType(args[1])
	} else {
		vp.SetConfigName("config")
		vp.SetConfigType("json")
	}
	vp.AddConfigPath("../")
	for _, path := range additionalPath {
		vp.AddConfigPath(path)
	}

	err := vp.ReadInConfig()
	if err != nil {
		log.Fatal("Cannot read the config file: ", err.Error())
		return Config{}
	}

	err = vp.Unmarshal(&config)
	if err != nil {
		log.Fatal("Cannot unmarshal the config file: ", err.Error())
		return Config{}
	}

	return config
}
